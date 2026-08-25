import {
  angleDegrees,
  cross,
  distance as coordinateDistance,
  subtract,
  vectorAngleDegrees,
} from "../geometry.js";
import type {
  AtomHybridization,
  ChemicalComponent,
  InteractionGeometryValue,
  InteractionMetadata,
  InteractionRule,
} from "../types.js";
import { SingleAngleRule } from "./angle.js";
import { InvertedRule } from "./inverted.js";
import { atomAt, interactionMetadata } from "./metadata.js";
import {
  elementSymbol,
  PROLIF_CSD_VDW_RADII,
} from "./vdw.js";

export const PROLIF_2_2_HBOND_ACCEPTOR_SMARTS =
  "[$([N&!$([NX3]-*=[O,N,P,S])&!$([ND2v3^2+0](-[H])-[CD4v4H1^3]-[CD2^2+0]=O)&!$([NX3]-[a])&!$([Nv4+1])&!$(N=C(-[C,N])-N)]),$([n+0&!X3&!$([n&r5]:[n+&r5])]),$([O&!$([OX2](C)C=O)&!$(O(~a)~a)&!$(O=N-*)&!$([O-]-N=O)]),$([o+0]),$([F&$(F-[#6])&!$(F-[#6][F,Cl,Br,I])])]";

export const PROLIF_2_2_HBOND_DONOR_SMARTS =
  "[$([O,S,#7;+0]),$([Nv4+1]),$([n+]c[nH])]-[H]";

export const PROLIF_2_2_IMPLICIT_HBOND_ACCEPTOR_SMARTS =
  PROLIF_2_2_HBOND_ACCEPTOR_SMARTS;

export const PROLIF_2_2_IMPLICIT_HBOND_DONOR_SMARTS =
  "[$([O,S,#7;+0&h,+0&H,+0&H2,+0H3]),$([N;v4&+1&h,v4&+1&H,v4&+1&H2,v4&+1&H3,v4&+1&H4]),$([n+]c[nh]),$([n+]c[nH])]";

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

const IDEAL_ATOM_ANGLES: Readonly<Record<AtomHybridization, number | undefined>> = {
  SP: 180,
  SP2: 120,
  SP3: 109.5,
  OTHER: undefined,
};

interface AtomGeometry {
  readonly idealAngle: number;
  readonly atomAngles: readonly number[];
  readonly deviation: number;
  readonly planeAngle?: number;
}

export interface ImplicitHBAcceptorOptions {
  readonly acceptor?: string;
  readonly donor?: string;
  readonly distance?: number;
  readonly toleranceDeviationDonorAtomAngle?: number;
  readonly toleranceDeviationDonorPlaneAngle?: number;
  readonly vinaPotentialMaximum?: number;
  readonly vinaPotentialMinimum?: number;
  readonly ignoreGeometryChecks?: boolean;
}

function heavyNeighbors(
  component: ChemicalComponent,
  atomIndex: number,
): readonly number[] {
  return component
    .neighbors(atomIndex)
    .filter((index) => atomAt(component, index, "ligand").atomicNumber !== 1);
}

function atomGeometry(
  component: ChemicalComponent,
  atomIndex: number,
  remoteComponent: ChemicalComponent,
  remoteAtomIndex: number,
): AtomGeometry {
  const center = atomAt(component, atomIndex, "ligand");
  const remote = atomAt(remoteComponent, remoteAtomIndex, "protein");
  const neighbors = [...heavyNeighbors(component, atomIndex)];
  if (neighbors.length === 0) {
    throw new RangeError(
      `No nearby heavy atoms found for atom index ${atomIndex}`,
    );
  }

  const hybridization = component.hybridization(atomIndex);
  const idealAngle = IDEAL_ATOM_ANGLES[hybridization];
  if (idealAngle === undefined) {
    throw new RangeError(
      `Implicit H-bond geometry does not support ${hybridization} hybridization`,
    );
  }
  const atomAngles = neighbors.map((neighborIndex) =>
    angleDegrees(
      atomAt(component, neighborIndex, "ligand").position,
      center.position,
      remote.position,
    ),
  );
  const deviation = Math.min(
    ...atomAngles.map((angle) => Math.abs(angle - idealAngle)),
  );

  if (hybridization !== "SP2") {
    return { idealAngle, atomAngles, deviation };
  }

  if (neighbors.length < 2) {
    for (const neighborIndex of [...neighbors]) {
      for (const extendedIndex of heavyNeighbors(component, neighborIndex)) {
        if (extendedIndex !== atomIndex) neighbors.push(extendedIndex);
      }
    }
  }
  const first = neighbors[0];
  const second = neighbors[1];
  if (first === undefined || second === undefined) {
    throw new RangeError(
      `Cannot construct an SP2 plane for atom index ${atomIndex}`,
    );
  }
  const normal = cross(
    subtract(atomAt(component, first, "ligand").position, center.position),
    subtract(atomAt(component, second, "ligand").position, center.position),
  );
  const planeAngle = Math.abs(
    vectorAngleDegrees(normal, subtract(remote.position, center.position)) - 90,
  );
  return { idealAngle, atomAngles, deviation, planeAngle };
}

/** Browser port of ProLIF 2.2.1's heavy-atom implicit H-bond rule. */
export class ImplicitHBAcceptor implements InteractionRule {
  readonly name = "ImplicitHBAcceptor";
  readonly acceptorSmarts: string;
  readonly donorSmarts: string;
  readonly cutoff: number;
  readonly toleranceDeviationDonorAtomAngle: number;
  readonly toleranceDeviationDonorPlaneAngle: number;
  readonly vinaPotentialMaximum: number;
  readonly vinaPotentialMinimum: number;
  readonly ignoreGeometryChecks: boolean;

  constructor(options: ImplicitHBAcceptorOptions = {}) {
    this.acceptorSmarts =
      options.acceptor ?? PROLIF_2_2_IMPLICIT_HBOND_ACCEPTOR_SMARTS;
    this.donorSmarts = options.donor ?? PROLIF_2_2_IMPLICIT_HBOND_DONOR_SMARTS;
    this.cutoff = options.distance ?? 3.5;
    this.toleranceDeviationDonorAtomAngle =
      options.toleranceDeviationDonorAtomAngle ?? 25;
    this.toleranceDeviationDonorPlaneAngle =
      options.toleranceDeviationDonorPlaneAngle ?? 30;
    this.vinaPotentialMaximum = options.vinaPotentialMaximum ?? -0.425;
    this.vinaPotentialMinimum = options.vinaPotentialMinimum ?? 0.565;
    this.ignoreGeometryChecks = options.ignoreGeometryChecks ?? false;
  }

  detect(
    ligand: ChemicalComponent,
    protein: ChemicalComponent,
  ): readonly InteractionMetadata[] {
    const output: InteractionMetadata[] = [];
    for (const ligandMatch of ligand.findMatches(this.acceptorSmarts)) {
      const ligandIndex = ligandMatch[0];
      if (ligandIndex === undefined) continue;
      const ligandAtom = atomAt(ligand, ligandIndex, "ligand");
      for (const proteinMatch of protein.findMatches(this.donorSmarts)) {
        const proteinIndex = proteinMatch[0];
        if (proteinIndex === undefined) continue;
        const proteinAtom = atomAt(protein, proteinIndex, "protein");
        const separation = coordinateDistance(
          ligandAtom.position,
          proteinAtom.position,
        );
        if (separation > this.cutoff) continue;

        const geometry: Record<string, InteractionGeometryValue> = {};
        if (!this.ignoreGeometryChecks) {
          const acceptor = atomGeometry(
            ligand,
            ligandIndex,
            protein,
            proteinIndex,
          );
          if (acceptor.deviation > 90) continue;
          geometry.ideal_acceptor_angle = acceptor.idealAngle;
          geometry.acceptor_atom_angles = acceptor.atomAngles;
          geometry.acceptor_atom_angle_deviation = acceptor.deviation;
          if (acceptor.planeAngle !== undefined) {
            if (acceptor.planeAngle > 90) continue;
            geometry.acceptor_plane_angle = acceptor.planeAngle;
          }

          const donor = atomGeometry(
            protein,
            proteinIndex,
            ligand,
            ligandIndex,
          );
          if (donor.deviation > this.toleranceDeviationDonorAtomAngle) continue;
          geometry.ideal_donor_angle = donor.idealAngle;
          geometry.donor_atom_angles = donor.atomAngles;
          geometry.donor_atom_angle_deviation = donor.deviation;
          if (donor.planeAngle !== undefined) {
            if (donor.planeAngle > this.toleranceDeviationDonorPlaneAngle) continue;
            geometry.donor_plane_angle = donor.planeAngle;
          }
        }

        const ligandRadius =
          PROLIF_CSD_VDW_RADII[elementSymbol(ligandAtom.atomicNumber)];
        const proteinRadius =
          PROLIF_CSD_VDW_RADII[elementSymbol(proteinAtom.atomicNumber)];
        if (ligandRadius === undefined || proteinRadius === undefined) {
          throw new RangeError("Missing CSD VdW radius for implicit H-bond atoms");
        }
        const distanceDifference = separation - ligandRadius - proteinRadius;
        geometry.vina_hbond_potential =
          distanceDifference <= this.vinaPotentialMaximum
            ? 1
            : distanceDifference >= this.vinaPotentialMinimum
              ? 0
              : (distanceDifference - this.vinaPotentialMinimum) /
                (this.vinaPotentialMaximum - this.vinaPotentialMinimum);

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
    return output;
  }
}

export class ImplicitHBDonor extends InvertedRule {
  constructor(options: ImplicitHBAcceptorOptions = {}) {
    super("ImplicitHBDonor", new ImplicitHBAcceptor(options));
  }
}
