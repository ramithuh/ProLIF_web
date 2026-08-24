import {
  betweenInclusive,
  centroid,
  distance,
  foldRingAngle,
  ringNormal,
  subtract,
  vectorAngleDegrees,
} from "../geometry.js";
import type {
  ChemicalComponent,
  InteractionMetadata,
  InteractionRule,
  Vec3,
} from "../types.js";
import { InvertedRule } from "./inverted.js";
import { PROLIF_2_2_CATION_SMARTS } from "./ionic.js";
import { atomAt, interactionMetadata } from "./metadata.js";

export const PROLIF_AROMATIC_RINGS = [
  "[a;r6]1:[a;r6]:[a;r6]:[a;r6]:[a;r6]:[a;r6]:1",
  "[a;r5]1:[a;r5]:[a;r5]:[a;r5]:[a;r5]:1",
] as const;

function coordinates(
  component: ChemicalComponent,
  indices: readonly number[],
  role: "ligand" | "protein",
): readonly Vec3[] {
  return indices.map((index) => atomAt(component, index, role).position);
}

export class CationPi implements InteractionRule {
  readonly name = "CationPi";

  constructor(
    readonly distanceCutoff = 4.5,
    readonly angle: readonly [number, number] = [0, 30],
  ) {}

  detect(
    ligand: ChemicalComponent,
    protein: ChemicalComponent,
  ): readonly InteractionMetadata[] {
    const output: InteractionMetadata[] = [];
    const cationMatches = ligand.findMatches(PROLIF_2_2_CATION_SMARTS);
    for (const ringSmarts of PROLIF_AROMATIC_RINGS) {
      for (const cationMatch of cationMatches) {
        const cationIndex = cationMatch[0];
        if (cationIndex === undefined) throw new RangeError("Cation match is empty");
        const cation = atomAt(ligand, cationIndex, "ligand");
        for (const ringMatch of protein.findMatches(ringSmarts)) {
          const ringCoordinates = coordinates(protein, ringMatch, "protein");
          const center = centroid(ringCoordinates);
          const separation = distance(cation.position, center);
          if (separation > this.distanceCutoff) continue;
          const normal = ringNormal(ringCoordinates);
          const measured = vectorAngleDegrees(
            normal,
            subtract(cation.position, center),
          );
          if (!betweenInclusive(foldRingAngle(measured), this.angle)) continue;
          output.push(
            interactionMetadata(
              this.name,
              ligand,
              protein,
              cationMatch,
              ringMatch,
              separation,
              { angle: measured },
            ),
          );
        }
      }
    }
    return output;
  }
}

export class PiCation extends InvertedRule {
  constructor(
    distanceCutoff = 4.5,
    angle: readonly [number, number] = [0, 30],
  ) {
    super("PiCation", new CationPi(distanceCutoff, angle));
  }
}
