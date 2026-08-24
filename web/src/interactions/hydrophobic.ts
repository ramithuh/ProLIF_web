import { DistanceRule } from "./distance.js";

/** ProLIF 2.2.x `Hydrophobic` SMARTS. */
export const PROLIF_2_2_HYDROPHOBIC_SMARTS =
  "[c,s,Br,I,S&H0&v2,$([C;X2H0R0,X3H1R0,X4H2R0,X3H0,X4H1,X4H0]);+0;!$([#6]~[#7,#8,#9])]";

export class Hydrophobic extends DistanceRule {
  constructor(distance = 4.5, smarts = PROLIF_2_2_HYDROPHOBIC_SMARTS) {
    super({
      name: "Hydrophobic",
      ligandSmarts: smarts,
      proteinSmarts: smarts,
      distance,
    });
  }
}
