import type {
  ChemicalComponent,
  InteractionMetadata,
  InteractionRule,
} from "./types.js";
import { CationPi, PiCation } from "./interactions/cation-pi.js";
import { HBAcceptor, HBDonor } from "./interactions/hbond.js";
import { Hydrophobic } from "./interactions/hydrophobic.js";
import { Anionic, Cationic } from "./interactions/ionic.js";
import { PiStacking } from "./interactions/pi-stacking.js";
import { VdWContact } from "./interactions/vdw.js";

export const DEFAULT_RULES = [
  new Hydrophobic(),
  new HBAcceptor(),
  new HBDonor(),
  new Cationic(),
  new Anionic(),
  new CationPi(),
  new PiCation(),
  new PiStacking(),
  new VdWContact(),
] as const satisfies readonly InteractionRule[];

/** Run a selected set of parity-tested rules over one ligand/protein pair. */
export class Fingerprint {
  constructor(readonly rules: readonly InteractionRule[] = DEFAULT_RULES) {}

  detect(
    ligand: ChemicalComponent,
    protein: ChemicalComponent,
  ): readonly InteractionMetadata[] {
    return this.rules.flatMap((rule) => rule.detect(ligand, protein));
  }
}

export { VdWContact };
