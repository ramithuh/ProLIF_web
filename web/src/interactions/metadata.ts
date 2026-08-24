import type {
  AtomMatch,
  AtomRecord,
  ChemicalComponent,
  InteractionMetadata,
} from "../types.js";

export function atomAt(
  component: ChemicalComponent,
  index: number,
  role: "ligand" | "protein",
): AtomRecord {
  const atom = component.atoms[index];
  if (atom === undefined) {
    throw new RangeError(`${role} match references missing atom index ${index}`);
  }
  return atom;
}

function mapped(
  component: ChemicalComponent,
  match: AtomMatch,
  role: "ligand" | "protein",
): readonly number[] {
  return match.map((index) => atomAt(component, index, role).parentIndex);
}

export function interactionMetadata(
  name: string,
  ligand: ChemicalComponent,
  protein: ChemicalComponent,
  ligandMatch: AtomMatch,
  proteinMatch: AtomMatch,
  distance: number,
  geometry?: Readonly<Record<string, number>>,
): InteractionMetadata {
  const common = {
    interaction: name,
    indices: { ligand: [...ligandMatch], protein: [...proteinMatch] },
    parentIndices: {
      ligand: mapped(ligand, ligandMatch, "ligand"),
      protein: mapped(protein, proteinMatch, "protein"),
    },
    distance,
  };
  return geometry === undefined ? common : { ...common, geometry };
}

export function invertedMetadata(
  name: string,
  metadata: InteractionMetadata,
): InteractionMetadata {
  return {
    ...metadata,
    interaction: name,
    indices: {
      ligand: metadata.indices.protein,
      protein: metadata.indices.ligand,
    },
    parentIndices: {
      ligand: metadata.parentIndices.protein,
      protein: metadata.parentIndices.ligand,
    },
  };
}
