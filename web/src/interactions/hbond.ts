import { SingleAngleRule } from "./angle.js";
import { InvertedRule } from "./inverted.js";

export const PROLIF_2_2_HBOND_ACCEPTOR_SMARTS =
  "[$([N&!$([NX3]-*=[O,N,P,S])&!$([ND2v3^2+0](-[H])-[CD4v4H1^3]-[CD2^2+0]=O)&!$([NX3]-[a])&!$([Nv4+1])&!$(N=C(-[C,N])-N)]),$([n+0&!X3&!$([n&r5]:[n+&r5])]),$([O&!$([OX2](C)C=O)&!$(O(~a)~a)&!$(O=N-*)&!$([O-]-N=O)]),$([o+0]),$([F&$(F-[#6])&!$(F-[#6][F,Cl,Br,I])])]";

export const PROLIF_2_2_HBOND_DONOR_SMARTS =
  "[$([O,S,#7;+0]),$([Nv4+1]),$([n+]c[nH])]-[H]";

export class HBAcceptor extends SingleAngleRule {
  constructor(
    distance = 3.5,
    angle: readonly [number, number] = [130, 180],
  ) {
    super({
      name: "HBAcceptor",
      ligandSmarts: PROLIF_2_2_HBOND_ACCEPTOR_SMARTS,
      proteinSmarts: PROLIF_2_2_HBOND_DONOR_SMARTS,
      distance,
      angle,
      distanceAtom: "P1",
      angleMetadataKey: "DHA_angle",
    });
  }
}

export class HBDonor extends InvertedRule {
  constructor(
    distance = 3.5,
    angle: readonly [number, number] = [130, 180],
  ) {
    super("HBDonor", new HBAcceptor(distance, angle));
  }
}
