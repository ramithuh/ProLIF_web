import type {
  AtomMatch,
  AtomHybridization,
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
    private readonly adjacency: Readonly<Record<number, readonly number[]>> = {},
    private readonly hybridizations: Readonly<Record<number, AtomHybridization>> = {},
  ) {}

  findMatches(smarts: string): readonly AtomMatch[] {
    return this.matches[smarts] ?? [];
  }

  neighbors(atomIndex: number): readonly number[] {
    return this.adjacency[atomIndex] ?? [];
  }

  hybridization(atomIndex: number): AtomHybridization {
    return this.hybridizations[atomIndex] ?? "OTHER";
  }
}
