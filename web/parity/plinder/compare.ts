export interface Event {
  interaction: string;
  indices: { ligand: number[]; protein: number[] };
  parentIndices: { ligand: number[]; protein: number[] };
  distance: number;
  geometry?: Record<string, unknown>;
}
const key = (e: Event) => JSON.stringify([e.interaction, e.indices.ligand, e.indices.protein,
  e.parentIndices.ligand, e.parentIndices.protein]);

function differences(a: unknown, b: unknown, path = '', tolerance = 1e-5): unknown[] {
  if (typeof a === 'number' && typeof b === 'number') {
    return Number.isFinite(a) && Number.isFinite(b) && Math.abs(a-b) <= tolerance
      ? [] : [{ path, native: a, web: b }];
  }
  if (a !== null && b !== null && typeof a === 'object' && typeof b === 'object') {
    if (Array.isArray(a) !== Array.isArray(b)) return [{ path, native: a, web: b }];
    const left = a as Record<string, unknown>, right = b as Record<string, unknown>;
    const keys = [...new Set([...Object.keys(left), ...Object.keys(right)])].sort();
    return keys.flatMap(k => differences(left[k], right[k], `${path}.${k}`, tolerance));
  }
  return a === b ? [] : [{ path, native: a ?? null, web: b ?? null }];
}

export function compareEvents(native: readonly Event[], web: readonly Event[], tolerance = 1e-5) {
  const remaining = [...web], nativeOnly: Event[] = [], geometryMismatches: unknown[] = [];
  let shared = 0;
  const diff = (a: Event, b: Event) => differences(
    { distance: a.distance, geometry: a.geometry ?? {} },
    { distance: b.distance, geometry: b.geometry ?? {} }, '', tolerance);
  for (const event of native) {
    const candidates = remaining.map((e,i) => ({e,i})).filter(({e}) => key(e) === key(event));
    const selected = candidates.find(({e}) => !diff(event,e).length) ?? candidates[0];
    if (!selected) { nativeOnly.push(event); continue; }
    shared++;
    remaining.splice(selected.i,1);
    const delta = diff(event,selected.e);
    if (delta.length) geometryMismatches.push({ event, differences: delta });
  }
  return { exact: !nativeOnly.length && !remaining.length && !geometryMismatches.length,
    nativeCount: native.length, webCount: web.length, shared,
    nativeOnly, webOnly: remaining, geometryMismatches };
}
