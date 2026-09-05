import type { JSMol, RDKitModule } from "@rdkit/rdkit";

import type {
  AtomHybridization,
  AtomMatch,
  AtomRecord,
  ChemicalComponent,
} from "../types.js";

interface RdkitMatch {
  readonly atoms?: unknown;
}

function parseMatches(
  payload: string,
  atomCount: number,
  smarts: string,
): readonly AtomMatch[] {
  let decoded: unknown;
  try {
    decoded = JSON.parse(payload);
  } catch (error) {
    throw new Error(`RDKit returned invalid match JSON for SMARTS ${smarts}`, {
      cause: error,
    });
  }
  // RDKit MinimalLib serializes an empty match set as `{}` rather than `[]`.
  if (
    typeof decoded === "object" &&
    decoded !== null &&
    !Array.isArray(decoded) &&
    Object.keys(decoded).length === 0
  ) {
    return [];
  }
  if (!Array.isArray(decoded)) {
    throw new TypeError(`RDKit match result for SMARTS ${smarts} is not an array`);
  }

  return decoded.map((rawMatch, matchIndex) => {
    const atoms = (rawMatch as RdkitMatch).atoms;
    if (!Array.isArray(atoms) || atoms.length === 0) {
      throw new TypeError(
        `RDKit match ${matchIndex} for SMARTS ${smarts} has no atoms`,
      );
    }
    return atoms.map((atomIndex) => {
      if (
        !Number.isInteger(atomIndex) ||
        (atomIndex as number) < 0 ||
        (atomIndex as number) >= atomCount
      ) {
        throw new RangeError(
          `RDKit match ${matchIndex} for SMARTS ${smarts} contains invalid atom index ${String(atomIndex)}`,
        );
      }
      return atomIndex as number;
    });
  });
}

/**
 * Adapts an RDKit.js molecule to the chemistry interface used by interaction
 * rules. Coordinates deliberately remain external: Mol*, PDB, and mmCIF
 * loaders can supply authoritative 3D positions without asking RDKit to alter
 * atom order or generate a conformer.
 */
export class RdkitChemicalComponent implements ChemicalComponent {
  readonly atoms: readonly AtomRecord[];
  private readonly matchCache = new Map<string, readonly AtomMatch[]>();
  private neighborCache: readonly (readonly number[])[] | undefined;
  private hybridizationCache: readonly AtomHybridization[] | undefined;
  private disposed = false;

  constructor(
    private readonly rdkit: RDKitModule,
    private readonly molecule: JSMol,
    atoms: readonly AtomRecord[],
    private readonly ownsMolecule = false,
    readonly residueName?: string,
  ) {
    atoms.forEach((atom, arrayIndex) => {
      if (atom.index !== arrayIndex) {
        throw new RangeError(
          `Atom array must follow RDKit ordering: position ${arrayIndex} has local index ${atom.index}`,
        );
      }
    });
    this.atoms = atoms;
  }

  findMatches(smarts: string): readonly AtomMatch[] {
    if (this.disposed) {
      throw new Error("Cannot query a disposed RDKit chemical component");
    }
    const cached = this.matchCache.get(smarts);
    if (cached !== undefined) {
      return cached;
    }

    const query = this.rdkit.get_qmol(smarts);
    if (query === null) {
      throw new SyntaxError(`RDKit could not parse SMARTS: ${smarts}`);
    }
    try {
      const matches = parseMatches(
        this.molecule.get_substruct_matches(query),
        this.atoms.length,
        smarts,
      );
      this.matchCache.set(smarts, matches);
      return matches;
    } finally {
      query.delete();
    }
  }

  neighbors(atomIndex: number): readonly number[] {
    this.assertAtomIndex(atomIndex);
    if (this.neighborCache === undefined) {
      const adjacency = Array.from(
        { length: this.atoms.length },
        () => [] as number[],
      );
      // GetNeighbors follows bond insertion order, not SMARTS match order.
      // Preserve that order: implicit H-bond plane construction uses it.
      const json = JSON.parse(this.molecule.get_json()) as {
        molecules: { bonds: { atoms: number[] }[] }[];
      };
      for (const bond of json.molecules[0]!.bonds) {
        const first = bond.atoms[0];
        const second = bond.atoms[1];
        if (first === undefined || second === undefined) {
          throw new RangeError("RDKit bond query returned an incomplete match");
        }
        adjacency[first]?.push(second);
        adjacency[second]?.push(first);
      }
      this.neighborCache = adjacency.map((indices) => [...indices]);
    }
    return this.neighborCache[atomIndex] ?? [];
  }

  hybridization(atomIndex: number): AtomHybridization {
    this.assertAtomIndex(atomIndex);
    if (this.hybridizationCache === undefined) {
      const values = Array.from<AtomHybridization>({ length: this.atoms.length }).fill(
        "OTHER",
      );
      for (const [smarts, hybridization] of [
        ["[*^1]", "SP"],
        ["[*^2]", "SP2"],
        ["[*^3]", "SP3"],
      ] as const) {
        for (const match of this.findMatches(smarts)) {
          const index = match[0];
          if (index !== undefined) values[index] = hybridization;
        }
      }
      this.hybridizationCache = values;
    }
    return this.hybridizationCache[atomIndex] ?? "OTHER";
  }

  private assertAtomIndex(atomIndex: number): void {
    if (!Number.isInteger(atomIndex) || atomIndex < 0 || atomIndex >= this.atoms.length) {
      throw new RangeError(`Invalid atom index ${atomIndex}`);
    }
  }

  /** Release an owned WebAssembly molecule. Safe to call more than once. */
  dispose(): void {
    if (!this.disposed && this.ownsMolecule) {
      this.molecule.delete();
    }
    this.disposed = true;
    this.matchCache.clear();
    this.neighborCache = undefined;
    this.hybridizationCache = undefined;
  }
}

/**
 * Parse a molecule with RDKit.js and return an owned component. The atom array
 * must use exactly the same atom ordering as the serialized molecule.
 */
export function componentFromRdkitInput(
  rdkit: RDKitModule,
  moleculeInput: string,
  atoms: readonly AtomRecord[],
  details?: Readonly<Record<string, unknown>>,
  context?: { readonly residueName?: string },
): RdkitChemicalComponent {
  // RDKit.js declares the second argument optional, but its Embind boundary
  // rejects an explicitly passed `undefined`. Invoke the one-argument form.
  const molecule =
    details === undefined
      ? rdkit.get_mol(moleculeInput)
      : rdkit.get_mol(moleculeInput, JSON.stringify(details));
  if (molecule === null) {
    throw new SyntaxError("RDKit could not parse the supplied molecule");
  }
  try {
    return new RdkitChemicalComponent(rdkit, molecule, atoms, true, context?.residueName);
  } catch (error) {
    molecule.delete();
    throw error;
  }
}
