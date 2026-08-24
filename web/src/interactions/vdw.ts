import { distance } from "../geometry.js";
import type {
  ChemicalComponent,
  InteractionMetadata,
  InteractionRule,
} from "../types.js";
import { interactionMetadata } from "./metadata.js";

const ELEMENTS = [
  "",
  "H",
  "He",
  "Li",
  "Be",
  "B",
  "C",
  "N",
  "O",
  "F",
  "Ne",
  "Na",
  "Mg",
  "Al",
  "Si",
  "P",
  "S",
  "Cl",
  "Ar",
  "K",
  "Ca",
  "Sc",
  "Ti",
  "V",
  "Cr",
  "Mn",
  "Fe",
  "Co",
  "Ni",
  "Cu",
  "Zn",
  "Ga",
  "Ge",
  "As",
  "Se",
  "Br",
  "Kr",
  "Rb",
  "Sr",
  "Y",
  "Zr",
  "Nb",
  "Mo",
  "Tc",
  "Ru",
  "Rh",
  "Pd",
  "Ag",
  "Cd",
  "In",
  "Sn",
  "Sb",
  "Te",
  "I",
  "Xe",
  "Cs",
  "Ba",
  "La",
  "Ce",
  "Pr",
  "Nd",
  "Pm",
  "Sm",
  "Eu",
  "Gd",
  "Tb",
  "Dy",
  "Ho",
  "Er",
  "Tm",
  "Yb",
  "Lu",
  "Hf",
  "Ta",
  "W",
  "Re",
  "Os",
  "Ir",
  "Pt",
  "Au",
  "Hg",
  "Tl",
  "Pb",
  "Bi",
  "Po",
  "At",
  "Rn",
  "Fr",
  "Ra",
  "Ac",
  "Th",
  "Pa",
  "U",
] as const;

/** ProLIF 2.2.1's default MDAnalysis radii preset. */
export const PROLIF_MDA_VDW_RADII: Readonly<Record<string, number>> = {
  H: 1.1,
  He: 1.4,
  Li: 1.82,
  Be: 1.53,
  B: 1.92,
  C: 1.7,
  N: 1.55,
  O: 1.52,
  F: 1.47,
  Ne: 1.54,
  Na: 2.27,
  Mg: 1.73,
  Al: 1.84,
  Si: 2.1,
  P: 1.8,
  S: 1.8,
  Cl: 1.75,
  Ar: 1.88,
  K: 2.75,
  Ca: 2.31,
  Ni: 1.63,
  Cu: 1.4,
  Zn: 1.39,
  Ga: 1.87,
  Ge: 2.11,
  Se: 1.9,
  Br: 1.85,
  Kr: 2.02,
  Sr: 2.49,
  Pd: 1.63,
  Ag: 1.72,
  Cd: 1.58,
  In: 1.93,
  Sn: 2.17,
  Sb: 2.06,
  Te: 2.06,
  I: 1.98,
  Xe: 2.16,
  Cs: 3.43,
  Ba: 2.68,
  Pt: 1.75,
  Au: 1.66,
  Tl: 1.96,
  Pb: 2.02,
  Bi: 2.07,
  Po: 1.97,
  At: 2.02,
  Rn: 2.2,
  Fr: 3.48,
  Ra: 2.83,
  U: 1.86,
};

function elementSymbol(atomicNumber: number): string {
  const symbol = ELEMENTS[atomicNumber];
  if (symbol === undefined || symbol === "") {
    throw new RangeError(`No element symbol for atomic number ${atomicNumber}`);
  }
  return symbol;
}

export class VdWContact implements InteractionRule {
  readonly name = "VdWContact";

  constructor(
    readonly tolerance = 0,
    readonly radii: Readonly<Record<string, number>> = PROLIF_MDA_VDW_RADII,
  ) {
    if (tolerance < 0) {
      throw new RangeError("VdW tolerance must be non-negative");
    }
  }

  detect(
    ligand: ChemicalComponent,
    protein: ChemicalComponent,
  ): readonly InteractionMetadata[] {
    const output: InteractionMetadata[] = [];
    for (const ligandAtom of ligand.atoms) {
      const ligandElement = elementSymbol(ligandAtom.atomicNumber);
      const ligandRadius = this.radii[ligandElement];
      if (ligandRadius === undefined) {
        throw new RangeError(`No VdW radius for ${ligandElement}`);
      }
      for (const proteinAtom of protein.atoms) {
        const proteinElement = elementSymbol(proteinAtom.atomicNumber);
        const proteinRadius = this.radii[proteinElement];
        if (proteinRadius === undefined) {
          throw new RangeError(`No VdW radius for ${proteinElement}`);
        }
        const separation = distance(ligandAtom.position, proteinAtom.position);
        if (separation > ligandRadius + proteinRadius + this.tolerance) continue;
        output.push(
          interactionMetadata(
            this.name,
            ligand,
            protein,
            [ligandAtom.index],
            [proteinAtom.index],
            separation,
          ),
        );
      }
    }
    return output;
  }
}
