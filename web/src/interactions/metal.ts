import { DistanceRule } from "./distance.js";
import { InvertedRule } from "./inverted.js";

export const PROLIF_2_2_METAL_SMARTS = "[Ca,Cd,Co,Cu,Fe,Mg,Mn,Ni,Zn]";
export const PROLIF_2_2_METAL_LIGAND_SMARTS =
  "[O,#7&!$([nX3])&!$([NX3]-*=[!#6])&!$([NX3]-[a])&!$([NX4]),-{1-};!+{1-}]";

export class MetalDonor extends DistanceRule {
  constructor(distance = 2.8) {
    super({
      name: "MetalDonor",
      ligandSmarts: PROLIF_2_2_METAL_SMARTS,
      proteinSmarts: PROLIF_2_2_METAL_LIGAND_SMARTS,
      distance,
    });
  }
}

export class MetalAcceptor extends InvertedRule {
  constructor(distance = 2.8) {
    super("MetalAcceptor", new MetalDonor(distance));
  }
}
