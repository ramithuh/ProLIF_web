import { createHash } from 'node:crypto';
import { readFileSync, writeFileSync, readdirSync, mkdirSync, copyFileSync } from 'node:fs';
import { execFileSync } from 'node:child_process';
import { expect, it } from 'vitest';
import type { RDKitLoader } from '@rdkit/rdkit';
import * as port from '../../src/index';
import { compareEvents, type Event } from './compare';

const root = process.env.PROLIF_PARITY_DIR;
if (!root) throw new Error('Set PROLIF_PARITY_DIR to the frozen native corpus');
const read = (path: string) => JSON.parse(readFileSync(path,'utf8'));
const write = (path: string, value: unknown) => writeFileSync(path,JSON.stringify(value,null,2)+'\n');
const hash = (path: string) => createHash('sha256').update(readFileSync(path)).digest('hex');
const rules: Record<string, port.InteractionRule> = {
  Hydrophobic: new port.Hydrophobic(), HBAcceptor: new port.HBAcceptor(), HBDonor: new port.HBDonor(),
  ImplicitHBAcceptor: new port.ImplicitHBAcceptor(), ImplicitHBDonor: new port.ImplicitHBDonor(),
  XBAcceptor: new port.XBAcceptor(), XBDonor: new port.XBDonor(), Cationic: new port.Cationic(), Anionic: new port.Anionic(),
  CationPi: new port.CationPi(), PiCation: new port.PiCation(), PiStacking: new port.PiStacking(),
  MetalDonor: new port.MetalDonor(), MetalAcceptor: new port.MetalAcceptor(), VdWContact: new port.VdWContact(),
  VdWContactCSD: new port.VdWContact(0,port.PROLIF_CSD_VDW_RADII),
};

it('matches native ProLIF on every frozen Plinder residue pair', async () => {
  const module = await import('@rdkit/rdkit');
  const rdkit = await (module.default as unknown as RDKitLoader)();
  const manifest = read(`${root}/manifest.json`), results: any[] = [];
  expect(manifest.cases.length,'empty corpora cannot establish parity').toBeGreaterThan(0);
  for (const caseInfo of manifest.cases) {
    const folder = `${root}/${caseInfo.pdb}`, native = read(`${folder}/native.json`);
    const row: any = { pdb: caseInfo.pdb, system: caseInfo.system, errors: [], failures: [], serializationChanges: [], unsupported: [], lanes: {} };
    results.push(row);
    if (native.error) { row.errors.push({stage:'native',error:native.error}); continue; }
    row.versions = native.versions;
    row.errors.push(...native.errors);
    expect(Object.keys(rules).sort()).toEqual([...native.rules].sort());
    for (const [file, expectedHash] of Object.entries(read(`${folder}/inputs.json`))) {
      expect(hash(`${folder}/${file}`),'frozen native input hash').toBe(expectedHash);
    }
    const components = new Map<string, ReturnType<typeof port.componentFromRdkitInput>>();
    try {
      for (const [id, data] of Object.entries(native.components) as [string, any][]) {
        try {
          expect(createHash('sha256').update(data.molblock).digest('hex')).toBe(data.graphSha256);
          components.set(id,port.componentFromRdkitInput(rdkit,data.molblock,data.atoms,{removeHs:false,sanitize:true},
            {residueName: data.residue.match(/^(.*?)\d+(?:\.|$)/)?.[1] ?? data.residue}));
        } catch(e) { row.errors.push({stage:'web-prepare',component:id,error:String(e)}); }
      }
      for (const pair of native.pairs) {
        const counts = row.lanes[pair.lane] ??= { pairs:0, checks:0, positiveChecks:0, exactChecks:0,
          nativeEvents:0, webEvents:0, nativeOnly:0, webOnly:0, geometryMismatches:0, byRule:{} };
        counts.pairs++;
        const ligand = components.get(pair.ligand), protein = components.get(pair.protein);
        if (!ligand || !protein) continue;
        for (const [name, rule] of Object.entries(rules)) {
          const expected = pair.rules[name], context = {lane:pair.lane,ligand:pair.ligand,protein:pair.protein,rule:name};
          if (expected.error) {
            // Native errors are not valid positive/negative oracle examples.
            // Still run Web to distinguish shared unsupported input from a
            // one-sided behavioral mismatch; don't bury these in match counts.
            try {
              const actual = rule.detect(ligand,protein);
              row.errors.push({...context,stage:'native-only-error',error:expected.error,webEvents:actual});
            } catch(e) { row.unsupported.push({...context,nativeError:expected.error,webError:String(e)}); }
            continue;
          }
          try {
            const actual = rule.detect(ligand,protein) as unknown as Event[];
            const comparison = compareEvents(expected.expected,actual);
            if (expected.beforeSerializationError) row.serializationChanges.push({...context,error:expected.beforeSerializationError});
            else {
              const serialization = compareEvents(expected.beforeSerialization,expected.expected,0.01);
              if (!serialization.exact) row.serializationChanges.push({...context,...serialization});
            }
            counts.checks++;
            counts.positiveChecks += Number(comparison.nativeCount > 0 || comparison.webCount > 0);
            counts.exactChecks += Number(comparison.exact);
            counts.nativeEvents += comparison.nativeCount; counts.webEvents += comparison.webCount;
            counts.nativeOnly += comparison.nativeOnly.length; counts.webOnly += comparison.webOnly.length;
            counts.geometryMismatches += comparison.geometryMismatches.length;
            const r = counts.byRule[name] ??= {checks:0,positiveChecks:0,failures:0,nativeEvents:0,webEvents:0};
            r.checks++; r.positiveChecks += Number(comparison.nativeCount>0 || comparison.webCount>0);
            r.failures += Number(!comparison.exact); r.nativeEvents += comparison.nativeCount; r.webEvents += comparison.webCount;
            if (!comparison.exact) row.failures.push({...context,...comparison});
          } catch(e) { row.errors.push({...context,stage:'web-rule',error:String(e)}); }
        }
      }
    } finally { for (const component of components.values()) component.dispose(); }
    row.exact = !row.errors.length && !row.failures.length && !row.unsupported.length && native.pairs.length>0;
    write(`${folder}/comparison.json`,row);
    console.log(caseInfo.pdb, row.exact?'MATCH':'FAIL',row.failures.length,'event/geometry mismatches;',row.errors.length,'errors;',row.unsupported.length,'shared unsupported');
  }
  const summary: any = {sampled:results.length,exactStructures:results.filter(r=>r.exact).length,
    structuresWithErrors:results.filter(r=>r.errors.length).length,
    mismatchingChecks:results.reduce((s,r)=>s+r.failures.length,0),
    sharedUnsupportedChecks:results.reduce((s,r)=>s+r.unsupported.length,0),
    errorChecks:results.reduce((s,r)=>s+r.errors.length,0),lanes:{}};
  for (const lane of ['implicit','explicit']) {
    const total: any = {};
    for (const key of ['pairs','checks','positiveChecks','exactChecks','nativeEvents','webEvents','nativeOnly','webOnly','geometryMismatches']) {
      total[key] = results.reduce((s,r)=>s+(r.lanes[lane]?.[key]??0),0);
    }
    total.byRule = {};
    for (const name of Object.keys(rules)) {
      total.byRule[name] = {};
      for (const key of ['checks','positiveChecks','failures','nativeEvents','webEvents']) {
        total.byRule[name][key] = results.reduce((s,r)=>s+(r.lanes[lane]?.byRule[name]?.[key]??0),0);
      }
    }
    summary.lanes[lane] = total;
  }
  const sourceHashes: Record<string,string> = {};
  const snapshot = (from: string,label:string) => {
    for (const item of readdirSync(from,{withFileTypes:true})) {
      const path = `${from}/${item.name}`, key = `${label}/${item.name}`;
      if (item.isDirectory()) { if(item.name!=='__pycache__') snapshot(path,key); }
      else if (/\.(ts|py|json)$/.test(item.name)) {
        mkdirSync(`${root}/source-snapshot/${label}`,{recursive:true});
        copyFileSync(path,`${root}/source-snapshot/${key}`); sourceHashes[key]=hash(path);
      }
    }
  };
  snapshot('src','prolif-web/src'); snapshot('parity/plinder','harness');
  write(`${root}/report.json`,{summary,manifest,provenance:{webRdkit:rdkit.version(),
    portCommit:execFileSync('git',['rev-parse','HEAD'],{encoding:'utf8'}).trim(),sourceHashes,
    wasmSha256:hash('node_modules/@rdkit/rdkit/dist/RDKit_minimal.wasm'),
    scope:'Shared serialized molecules, coordinates and hydrogen states; not Weaver CCD preparation, browser UI or full ProLIF high-level pipeline.'},results});
  console.log(JSON.stringify(summary,null,2));
  if (process.env.PROLIF_PARITY_REPORT_ONLY!=='1') {
    expect(results.filter(r=>!r.exact).map(r=>r.pdb),'native ProLIF conformance failures').toEqual([]);
  }
},180000);
