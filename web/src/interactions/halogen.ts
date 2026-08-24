import { DoubleAngleRule } from "./angle.js";
import { InvertedRule } from "./inverted.js";

export const PROLIF_2_2_HALOGEN_ACCEPTOR_SMARTS =
  "[#7,#8,P,S,Se,Te,a;!+{1-}]!#[*]";
export const PROLIF_2_2_HALOGEN_DONOR_SMARTS =
  "[#6,#7,Si,F,Cl,Br,I]-[Cl,Br,I,At]";

export class XBAcceptor extends DoubleAngleRule {
  constructor(
    distance = 3.5,
    AXDangle: readonly [number, number] = [130, 180],
    XARangle: readonly [number, number] = [80, 140],
  ) {
    super({
      name: "XBAcceptor",
      ligandSmarts: PROLIF_2_2_HALOGEN_ACCEPTOR_SMARTS,
      proteinSmarts: PROLIF_2_2_HALOGEN_DONOR_SMARTS,
      distance,
      firstAngle: AXDangle,
      secondAngle: XARangle,
      distanceAtoms: ["L1", "P2"],
      firstAngleMetadataKey: "AXD_angle",
      secondAngleMetadataKey: "XAR_angle",
    });
  }
}

export class XBDonor extends InvertedRule {
  constructor(
    distance = 3.5,
    AXDangle: readonly [number, number] = [130, 180],
    XARangle: readonly [number, number] = [80, 140],
  ) {
    super("XBDonor", new XBAcceptor(distance, AXDangle, XARangle));
  }
}
