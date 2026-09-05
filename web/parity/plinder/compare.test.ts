import { describe, expect, it } from 'vitest';
import { compareEvents, type Event } from './compare';

const event = (change: Partial<Event> = {}): Event => ({interaction:'ImplicitHBAcceptor',
  indices:{ligand:[1],protein:[2]},parentIndices:{ligand:[11],protein:[22]},
  distance:3,geometry:{angles:[110,120]},...change});

describe('ProLIF conformance comparator',()=>{
  it('matches reordered events without dropping duplicates',()=>{
    const a=event(), b=event({distance:4});
    expect(compareEvents([a,b],[b,a]).exact).toBe(true);
    expect(compareEvents([a,a],[a]).nativeOnly).toHaveLength(1);
    expect(compareEvents([a],[a,a]).webOnly).toHaveLength(1);
  });
  it.each([
    {interaction:'ImplicitHBDonor'}, {indices:{ligand:[3],protein:[2]}},
    {parentIndices:{ligand:[11],protein:[99]}},
  ])('does not merge different atom-event identities: %j',change=>{
    const diff=compareEvents([event()],[event(change)]);
    expect(diff.nativeOnly).toHaveLength(1); expect(diff.webOnly).toHaveLength(1);
  });
  it('distinguishes geometry-array ordering from detection differences',()=>{
    const diff=compareEvents([event()],[event({geometry:{angles:[120,110]}})]);
    expect(diff.shared).toBe(1); expect(diff.nativeOnly).toHaveLength(0);
    expect(diff.geometryMismatches).toHaveLength(1);
  });
  it('rejects missing geometry, extra geometry, and nonfinite numbers',()=>{
    for(const geometry of [{},{angles:[110,120],extra:1},{angles:[110,NaN]}]) {
      expect(compareEvents([event()],[event({geometry})]).exact).toBe(false);
    }
    expect(compareEvents([event()],[event({distance:Infinity})]).exact).toBe(false);
  });
  it('allows floating point noise but not a materially different distance',()=>{
    expect(compareEvents([event()],[event({distance:3.000001})]).exact).toBe(true);
    expect(compareEvents([event()],[event({distance:3.001})]).exact).toBe(false);
  });
  it('compares empty and omitted geometry equally',()=>{
    const a=event({geometry:{}}); const b=event(); delete b.geometry;
    expect(compareEvents([a],[b]).exact).toBe(true);
  });
});
