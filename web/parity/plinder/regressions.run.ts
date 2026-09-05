import { readFileSync } from 'node:fs';
import { beforeAll, expect, it } from 'vitest';
import type { RDKitLoader, RDKitModule } from '@rdkit/rdkit';
import { componentFromRdkitInput, ImplicitHBAcceptor, ImplicitHBDonor } from '../../src/index';
import { compareEvents, type Event } from './compare';

const fixture = JSON.parse(readFileSync(new URL('./fixtures/regressions.json',import.meta.url),'utf8'));
let rdkit: RDKitModule;
beforeAll(async()=>{const module=await import('@rdkit/rdkit');rdkit=await (module.default as unknown as RDKitLoader)();});
for(const c of fixture.cases) {
  it(`strict native ProLIF regression ${c.pdb}: ${c.kind}`,()=>{
    const context=(data: {residue:string})=>({residueName:data.residue.match(/^(.*?)\d+(?:\.|$)/)?.[1] ?? data.residue});
    const ligand=componentFromRdkitInput(rdkit,c.ligand.molblock,c.ligand.atoms,{removeHs:false,sanitize:true},context(c.ligand));
    const protein=componentFromRdkitInput(rdkit,c.protein.molblock,c.protein.atoms,{removeHs:false,sanitize:true},context(c.protein));
    try {
      const rule=c.rule==='ImplicitHBDonor'?new ImplicitHBDonor():new ImplicitHBAcceptor();
      const actual=rule.detect(ligand,protein) as unknown as Event[];
      expect(compareEvents(c.expected,actual)).toMatchObject({exact:true});
    } finally {ligand.dispose();protein.dispose();}
  });
}
