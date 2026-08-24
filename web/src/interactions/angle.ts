import { angleDegrees, betweenInclusive, distance } from "../geometry.js";
import type {
  ChemicalComponent,
  InteractionMetadata,
  InteractionRule,
} from "../types.js";
import { atomAt, interactionMetadata } from "./metadata.js";

export interface SingleAngleRuleOptions {
  readonly name: string;
  readonly ligandSmarts: string;
  readonly proteinSmarts: string;
  readonly distance: number;
  readonly angle: readonly [number, number];
  readonly distanceAtom: "P1" | "P2";
  readonly angleMetadataKey: string;
}

export class SingleAngleRule implements InteractionRule {
  readonly name: string;

  constructor(readonly options: SingleAngleRuleOptions) {
    this.name = options.name;
  }

  detect(
    ligand: ChemicalComponent,
    protein: ChemicalComponent,
  ): readonly InteractionMetadata[] {
    const output: InteractionMetadata[] = [];
    for (const ligandMatch of ligand.findMatches(this.options.ligandSmarts)) {
      const l1Index = ligandMatch[0];
      if (l1Index === undefined) throw new RangeError("Ligand match is empty");
      const l1 = atomAt(ligand, l1Index, "ligand");
      for (const proteinMatch of protein.findMatches(this.options.proteinSmarts)) {
        const p1Index = proteinMatch[0];
        const p2Index = proteinMatch[1];
        if (p1Index === undefined || p2Index === undefined) {
          throw new RangeError("Protein angle match requires two atoms");
        }
        const p1 = atomAt(protein, p1Index, "protein");
        const p2 = atomAt(protein, p2Index, "protein");
        const separation = distance(
          l1.position,
          this.options.distanceAtom === "P1" ? p1.position : p2.position,
        );
        if (separation > this.options.distance) continue;
        const angle = angleDegrees(p1.position, p2.position, l1.position);
        if (!betweenInclusive(angle, this.options.angle)) continue;
        output.push(
          interactionMetadata(
            this.name,
            ligand,
            protein,
            ligandMatch,
            proteinMatch,
            separation,
            { [this.options.angleMetadataKey]: angle },
          ),
        );
      }
    }
    return output;
  }
}

export interface DoubleAngleRuleOptions {
  readonly name: string;
  readonly ligandSmarts: string;
  readonly proteinSmarts: string;
  readonly distance: number;
  readonly firstAngle: readonly [number, number];
  readonly secondAngle: readonly [number, number];
  readonly distanceAtoms: readonly ["L1" | "L2", "P1" | "P2"];
  readonly firstAngleMetadataKey: string;
  readonly secondAngleMetadataKey: string;
}

export class DoubleAngleRule implements InteractionRule {
  readonly name: string;

  constructor(readonly options: DoubleAngleRuleOptions) {
    this.name = options.name;
  }

  detect(
    ligand: ChemicalComponent,
    protein: ChemicalComponent,
  ): readonly InteractionMetadata[] {
    const output: InteractionMetadata[] = [];
    for (const ligandMatch of ligand.findMatches(this.options.ligandSmarts)) {
      const l1Index = ligandMatch[0];
      const l2Index = ligandMatch[1];
      if (l1Index === undefined || l2Index === undefined) {
        throw new RangeError("Ligand angle match requires two atoms");
      }
      const l1 = atomAt(ligand, l1Index, "ligand");
      const l2 = atomAt(ligand, l2Index, "ligand");
      for (const proteinMatch of protein.findMatches(this.options.proteinSmarts)) {
        const p1Index = proteinMatch[0];
        const p2Index = proteinMatch[1];
        if (p1Index === undefined || p2Index === undefined) {
          throw new RangeError("Protein angle match requires two atoms");
        }
        const p1 = atomAt(protein, p1Index, "protein");
        const p2 = atomAt(protein, p2Index, "protein");
        const ligandDistanceAtom =
          this.options.distanceAtoms[0] === "L1" ? l1 : l2;
        const proteinDistanceAtom =
          this.options.distanceAtoms[1] === "P1" ? p1 : p2;
        const separation = distance(
          ligandDistanceAtom.position,
          proteinDistanceAtom.position,
        );
        if (separation > this.options.distance) continue;
        const firstAngle = angleDegrees(p1.position, p2.position, l1.position);
        if (!betweenInclusive(firstAngle, this.options.firstAngle)) continue;
        const secondAngle = angleDegrees(l2.position, l1.position, p2.position);
        if (!betweenInclusive(secondAngle, this.options.secondAngle)) continue;
        output.push(
          interactionMetadata(
            this.name,
            ligand,
            protein,
            ligandMatch,
            proteinMatch,
            separation,
            {
              [this.options.firstAngleMetadataKey]: firstAngle,
              [this.options.secondAngleMetadataKey]: secondAngle,
            },
          ),
        );
      }
    }
    return output;
  }
}
