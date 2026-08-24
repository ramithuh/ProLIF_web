import type {
  AtomMatch,
  AtomRecord,
  ChemicalComponent,
  Vec3,
} from "../src/index.js";

export function atom(
  index: number,
  parentIndex: number,
  position: Vec3,
  atomicNumber = 6,
): AtomRecord {
  return { index, parentIndex, atomicNumber, position };
}

export class FixtureComponent implements ChemicalComponent {
  constructor(
    readonly atoms: readonly AtomRecord[],
    private readonly matches: Readonly<Record<string, readonly AtomMatch[]>>,
  ) {}

  findMatches(smarts: string): readonly AtomMatch[] {
    return this.matches[smarts] ?? [];
  }
}
