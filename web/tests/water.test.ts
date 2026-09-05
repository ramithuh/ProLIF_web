import { beforeAll, expect, it } from 'vitest';
import type { RDKitLoader, RDKitModule } from '@rdkit/rdkit';
import { componentFromRdkitInput, ImplicitHBAcceptor, ImplicitHBDonor } from '../src/index';
let rdkit: RDKitModule;
beforeAll(async()=>{const module=await import('@rdkit/rdkit');rdkit=await (module.default as unknown as RDKitLoader)();});
it.each(['HOH','H20','WAT','SOL','TIP3','TP3','TIP'])('matches native water exclusion for %s',residueName=>{
  const make=(x:number)=>componentFromRdkitInput(rdkit,'O',[{index:0,parentIndex:0,atomicNumber:8,position:[x,0,0]}],
    {removeHs:false},{residueName});
  const a=make(0),b=make(2.8);
  try {
    expect(new ImplicitHBAcceptor().detect(a,b)).toEqual([]);
    expect(new ImplicitHBDonor().detect(a,b)).toEqual([]);
    expect(new ImplicitHBAcceptor({includeWater:true}).detect(a,b)).toHaveLength(1);
    expect(new ImplicitHBDonor({includeWater:true}).detect(a,b)).toHaveLength(1);
  } finally {a.dispose();b.dispose();}
});
