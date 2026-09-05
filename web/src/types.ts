/** A Cartesian coordinate in Ångströms. */
export type Vec3 = readonly [x: number, y: number, z: number];

/**
 * One atom in the local molecule presented to an interaction rule.
 *
 * `index` mirrors ProLIF's residue-local index. `parentIndex` identifies the
 * atom in the unsplit parent molecule, matching ProLIF's `parent_indices`.
 */
export interface AtomRecord {
  readonly index: number;
  readonly parentIndex: number;
  readonly atomicNumber: number;
  readonly position: Vec3;
}

/** Ordered atom indices returned by one SMARTS substructure match. */
export type AtomMatch = readonly number[];

/** Hybridization states used by ProLIF's implicit-hydrogen geometry checks. */
export type AtomHybridization = "SP" | "SP2" | "SP3" | "OTHER";

/** Minimum interface required from a chemistry backend such as RDKit.js. */
export interface ChemicalComponent {
  /** Native residue identity, required for residue-aware rules such as water exclusion. */
  readonly residueName?: string | undefined;
  readonly atoms: readonly AtomRecord[];
  findMatches(smarts: string): readonly AtomMatch[];
  /** Atom indices directly bonded to `atomIndex`. */
  neighbors(atomIndex: number): readonly number[];
  /** RDKit hybridization for one atom. */
  hybridization(atomIndex: number): AtomHybridization;
}

export interface InteractionAtomIndices {
  readonly ligand: readonly number[];
  readonly protein: readonly number[];
}

/**
 * Browser equivalent of ProLIF's interaction metadata.
 *
 * Rule-specific geometry is retained under `geometry` while the common
 * distance remains directly accessible for sorting and display. Implicit
 * hydrogen-bond rules also report the complete set of heavy-neighbor angles.
 */
export type InteractionGeometryValue = number | readonly number[];

export interface InteractionMetadata {
  readonly interaction: string;
  readonly indices: InteractionAtomIndices;
  readonly parentIndices: InteractionAtomIndices;
  readonly distance: number;
  readonly geometry?: Readonly<Record<string, InteractionGeometryValue>>;
}

export interface InteractionRule {
  readonly name: string;
  detect(
    ligand: ChemicalComponent,
    protein: ChemicalComponent,
  ): readonly InteractionMetadata[];
}
