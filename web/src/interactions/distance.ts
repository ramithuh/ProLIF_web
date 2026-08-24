import { distance } from "../geometry.js";
import type {
  AtomMatch,
  AtomRecord,
  ChemicalComponent,
  InteractionMetadata,
  InteractionRule,
} from "../types.js";

export interface DistanceRuleOptions {
  readonly name: string;
  readonly ligandSmarts: string;
  readonly proteinSmarts: string;
  readonly distance: number;
}

function firstAtom(
  component: ChemicalComponent,
  match: AtomMatch,
  role: "ligand" | "protein",
): AtomRecord {
  const atomIndex = match[0];
  if (atomIndex === undefined) {
    throw new RangeError(`${role} SMARTS match is empty`);
  }
  const atom = component.atoms[atomIndex];
  if (atom === undefined) {
    throw new RangeError(
      `${role} SMARTS match references missing atom index ${atomIndex}`,
    );
  }
  return atom;
}

function parentIndices(
  component: ChemicalComponent,
  match: AtomMatch,
  role: "ligand" | "protein",
): readonly number[] {
  return match.map((atomIndex) => {
    const atom = component.atoms[atomIndex];
    if (atom === undefined) {
      throw new RangeError(
        `${role} SMARTS match references missing atom index ${atomIndex}`,
      );
    }
    return atom.parentIndex;
  });
}

/** Browser equivalent of `prolif.interactions.base.Distance`. */
export class DistanceRule implements InteractionRule {
  readonly name: string;
  readonly ligandSmarts: string;
  readonly proteinSmarts: string;
  readonly cutoff: number;

  constructor(options: DistanceRuleOptions) {
    if (!Number.isFinite(options.distance) || options.distance < 0) {
      throw new RangeError("Distance cutoff must be a finite non-negative number");
    }
    this.name = options.name;
    this.ligandSmarts = options.ligandSmarts;
    this.proteinSmarts = options.proteinSmarts;
    this.cutoff = options.distance;
  }

  detect(
    ligand: ChemicalComponent,
    protein: ChemicalComponent,
  ): readonly InteractionMetadata[] {
    const ligandMatches = ligand.findMatches(this.ligandSmarts);
    const proteinMatches = protein.findMatches(this.proteinSmarts);
    const interactions: InteractionMetadata[] = [];

    for (const ligandMatch of ligandMatches) {
      const ligandAtom = firstAtom(ligand, ligandMatch, "ligand");
      for (const proteinMatch of proteinMatches) {
        const proteinAtom = firstAtom(protein, proteinMatch, "protein");
        const separation = distance(ligandAtom.position, proteinAtom.position);
        if (separation <= this.cutoff) {
          interactions.push({
            interaction: this.name,
            indices: {
              ligand: [...ligandMatch],
              protein: [...proteinMatch],
            },
            parentIndices: {
              ligand: parentIndices(ligand, ligandMatch, "ligand"),
              protein: parentIndices(protein, proteinMatch, "protein"),
            },
            distance: separation,
          });
        }
      }
    }
    return interactions;
  }
}
