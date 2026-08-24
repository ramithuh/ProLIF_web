import { DistanceRule } from "./distance.js";
import { InvertedRule } from "./inverted.js";

export const PROLIF_2_2_CATION_SMARTS =
  "[+{1-}!$(*~[*-{1-}]),$([NX3&!$([NX3]-O)]-[C]=[NX3+])]";
export const PROLIF_2_2_ANION_SMARTS =
  "[-{1-}!$(*~[*+{1-}]),$(O=[C,S,P]-[O-])]";

export class Cationic extends DistanceRule {
  constructor(distance = 4.5) {
    super({
      name: "Cationic",
      ligandSmarts: PROLIF_2_2_CATION_SMARTS,
      proteinSmarts: PROLIF_2_2_ANION_SMARTS,
      distance,
    });
  }
}

export class Anionic extends InvertedRule {
  constructor(distance = 4.5) {
    super("Anionic", new Cationic(distance));
  }
}
