import {
  betweenInclusive,
  centroid,
  distance,
  foldRingAngle,
  ringNormal,
  ringPlaneIntersection,
  subtract,
  vectorAngleDegrees,
} from "../geometry.js";
import type {
  ChemicalComponent,
  InteractionMetadata,
  InteractionRule,
  Vec3,
} from "../types.js";
import { PROLIF_AROMATIC_RINGS } from "./cation-pi.js";
import { atomAt, interactionMetadata } from "./metadata.js";

interface PiStackingGeometry {
  readonly name: string;
  readonly distance: number;
  readonly planeAngle: readonly [number, number];
  readonly normalToCentroidAngle: readonly [number, number];
  readonly intersect: boolean;
  readonly intersectRadius: number;
}

function coordinates(
  component: ChemicalComponent,
  indices: readonly number[],
  role: "ligand" | "protein",
): readonly Vec3[] {
  return indices.map((index) => atomAt(component, index, role).position);
}

class PiStackingGeometryRule implements InteractionRule {
  readonly name: string;

  constructor(readonly options: PiStackingGeometry) {
    this.name = options.name;
  }

  detect(
    ligand: ChemicalComponent,
    protein: ChemicalComponent,
  ): readonly InteractionMetadata[] {
    const output: InteractionMetadata[] = [];
    for (const proteinRingSmarts of PROLIF_AROMATIC_RINGS) {
      for (const ligandRingSmarts of PROLIF_AROMATIC_RINGS) {
        const proteinMatches = protein.findMatches(proteinRingSmarts);
        const ligandMatches = ligand.findMatches(ligandRingSmarts);
        for (const ligandMatch of ligandMatches) {
          const ligandCoordinates = coordinates(ligand, ligandMatch, "ligand");
          const ligandCentroid = centroid(ligandCoordinates);
          const ligandNormal = ringNormal(ligandCoordinates);
          for (const proteinMatch of proteinMatches) {
            const proteinCoordinates = coordinates(
              protein,
              proteinMatch,
              "protein",
            );
            const proteinCentroid = centroid(proteinCoordinates);
            const separation = distance(ligandCentroid, proteinCentroid);
            if (separation > this.options.distance) continue;
            const proteinNormal = ringNormal(proteinCoordinates);
            const planeAngle = vectorAngleDegrees(ligandNormal, proteinNormal);
            if (
              !betweenInclusive(
                foldRingAngle(planeAngle),
                this.options.planeAngle,
              )
            ) {
              continue;
            }

            const ligandNcc = vectorAngleDegrees(
              ligandNormal,
              subtract(proteinCentroid, ligandCentroid),
            );
            const proteinNcc = vectorAngleDegrees(
              proteinNormal,
              subtract(ligandCentroid, proteinCentroid),
            );
            const normalToCentroidAngle = betweenInclusive(
              foldRingAngle(ligandNcc),
              this.options.normalToCentroidAngle,
            )
              ? ligandNcc
              : betweenInclusive(
                    foldRingAngle(proteinNcc),
                    this.options.normalToCentroidAngle,
                  )
                ? proteinNcc
                : null;
            if (normalToCentroidAngle === null) continue;

            const geometry: Record<string, number> = {
              plane_angle: planeAngle,
              normal_to_centroid_angle: normalToCentroidAngle,
            };
            if (this.options.intersect) {
              const intersection = ringPlaneIntersection(
                ligandNormal,
                ligandCentroid,
                proteinNormal,
                proteinCentroid,
              );
              if (intersection === null) continue;
              const intersectionDistance = Math.min(
                distance(ligandCentroid, intersection),
                distance(proteinCentroid, intersection),
              );
              if (intersectionDistance > this.options.intersectRadius) continue;
              geometry.intersect_distance = intersectionDistance;
            }

            output.push(
              interactionMetadata(
                this.name,
                ligand,
                protein,
                ligandMatch,
                proteinMatch,
                separation,
                geometry,
              ),
            );
          }
        }
      }
    }
    return output;
  }
}

export class FaceToFace extends PiStackingGeometryRule {
  constructor(name = "FaceToFace") {
    super({
      name,
      distance: 5.5,
      planeAngle: [0, 35],
      normalToCentroidAngle: [0, 33],
      intersect: false,
      intersectRadius: 1.5,
    });
  }
}

export class EdgeToFace extends PiStackingGeometryRule {
  constructor(name = "EdgeToFace") {
    super({
      name,
      distance: 6.5,
      planeAngle: [50, 90],
      normalToCentroidAngle: [0, 30],
      intersect: true,
      intersectRadius: 1.5,
    });
  }
}

export class PiStacking implements InteractionRule {
  readonly name = "PiStacking";
  private readonly faceToFace = new FaceToFace(this.name);
  private readonly edgeToFace = new EdgeToFace(this.name);

  detect(
    ligand: ChemicalComponent,
    protein: ChemicalComponent,
  ): readonly InteractionMetadata[] {
    return [
      ...this.faceToFace.detect(ligand, protein),
      ...this.edgeToFace.detect(ligand, protein),
    ];
  }
}
