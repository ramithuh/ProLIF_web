import { readFileSync } from "node:fs";

import type { RDKitLoader, RDKitModule } from "@rdkit/rdkit";
import { beforeAll, describe, expect, it } from "vitest";

import {
  Anionic,
  Cationic,
  CationPi,
  componentFromRdkitInput,
  HBAcceptor,
  HBDonor,
  Hydrophobic,
  ImplicitHBAcceptor,
  ImplicitHBDonor,
  MetalAcceptor,
  MetalDonor,
  PiCation,
  PiStacking,
  PROLIF_2_2_ANION_SMARTS,
  PROLIF_2_2_CATION_SMARTS,
  PROLIF_2_2_HALOGEN_ACCEPTOR_SMARTS,
  PROLIF_2_2_HALOGEN_DONOR_SMARTS,
  PROLIF_2_2_HBOND_ACCEPTOR_SMARTS,
  PROLIF_2_2_HBOND_DONOR_SMARTS,
  PROLIF_2_2_HYDROPHOBIC_SMARTS,
  PROLIF_2_2_IMPLICIT_HBOND_ACCEPTOR_SMARTS,
  PROLIF_2_2_IMPLICIT_HBOND_DONOR_SMARTS,
  PROLIF_2_2_METAL_LIGAND_SMARTS,
  PROLIF_2_2_METAL_SMARTS,
  PROLIF_AROMATIC_RINGS,
  VdWContact,
  XBAcceptor,
  XBDonor,
  type AtomRecord,
  type InteractionMetadata,
  type InteractionRule,
} from "../src/index.js";

interface SerializedComponent {
  readonly molblock: string;
  readonly atoms: readonly AtomRecord[];
}

interface ParityCase {
  readonly name: string;
  readonly rule: string;
  readonly ligand: SerializedComponent;
  readonly protein: SerializedComponent;
  readonly expected: readonly InteractionMetadata[];
}

interface SmartsParityCase {
  readonly name: string;
  readonly patternKey: string;
  readonly pythonSmarts: readonly string[];
  readonly component: SerializedComponent;
  readonly expectedMatches: readonly (readonly number[])[];
}

interface ParityFixture {
  readonly generatedBy: string;
  readonly pythonRdkitVersion: string;
  readonly rules: readonly string[];
  readonly cases: readonly ParityCase[];
  readonly smartsCases: readonly SmartsParityCase[];
}

const RULES: Readonly<Record<string, () => InteractionRule>> = {
  Hydrophobic: () => new Hydrophobic(),
  HBAcceptor: () => new HBAcceptor(),
  HBDonor: () => new HBDonor(),
  ImplicitHBAcceptor: () => new ImplicitHBAcceptor(),
  ImplicitHBDonor: () => new ImplicitHBDonor(),
  XBAcceptor: () => new XBAcceptor(),
  XBDonor: () => new XBDonor(),
  Cationic: () => new Cationic(),
  Anionic: () => new Anionic(),
  CationPi: () => new CationPi(),
  PiCation: () => new PiCation(),
  PiStacking: () => new PiStacking(),
  MetalDonor: () => new MetalDonor(),
  MetalAcceptor: () => new MetalAcceptor(),
  VdWContact: () => new VdWContact(),
};

const SMARTS_PATTERNS: Readonly<Record<string, readonly string[]>> = {
  "Hydrophobic.lig_pattern": [PROLIF_2_2_HYDROPHOBIC_SMARTS],
  "HBAcceptor.lig_pattern": [PROLIF_2_2_HBOND_ACCEPTOR_SMARTS],
  "HBAcceptor.prot_pattern": [PROLIF_2_2_HBOND_DONOR_SMARTS],
  "ImplicitHBAcceptor.lig_pattern": [
    PROLIF_2_2_IMPLICIT_HBOND_ACCEPTOR_SMARTS,
  ],
  "ImplicitHBAcceptor.prot_pattern": [
    PROLIF_2_2_IMPLICIT_HBOND_DONOR_SMARTS,
  ],
  "XBAcceptor.lig_pattern": [PROLIF_2_2_HALOGEN_ACCEPTOR_SMARTS],
  "XBAcceptor.prot_pattern": [PROLIF_2_2_HALOGEN_DONOR_SMARTS],
  "Cationic.lig_pattern": [PROLIF_2_2_CATION_SMARTS],
  "Cationic.prot_pattern": [PROLIF_2_2_ANION_SMARTS],
  "CationPi.cation": [PROLIF_2_2_CATION_SMARTS],
  "CationPi.pi_ring": PROLIF_AROMATIC_RINGS,
  "MetalDonor.lig_pattern": [PROLIF_2_2_METAL_SMARTS],
  "MetalDonor.prot_pattern": [PROLIF_2_2_METAL_LIGAND_SMARTS],
};

async function loadRdkit(): Promise<RDKitModule> {
  const namespace: unknown = await import("@rdkit/rdkit");
  const loader = (namespace as { readonly default?: RDKitLoader }).default;
  if (typeof loader !== "function") {
    throw new TypeError("@rdkit/rdkit did not expose its module loader");
  }
  return loader();
}

const fixture = JSON.parse(
  readFileSync(
    new URL("./reference/prolif-2.2.1-core.json", import.meta.url),
    "utf8",
  ),
) as ParityFixture;

function eventIdentity(event: InteractionMetadata): string {
  return JSON.stringify({
    interaction: event.interaction,
    indices: event.indices,
    parentIndices: event.parentIndices,
  });
}

function multisetIntersectionSize(
  left: readonly InteractionMetadata[],
  right: readonly InteractionMetadata[],
): number {
  const available = new Map<string, number>();
  for (const event of right) {
    const key = eventIdentity(event);
    available.set(key, (available.get(key) ?? 0) + 1);
  }
  let matched = 0;
  for (const event of left) {
    const key = eventIdentity(event);
    const count = available.get(key) ?? 0;
    if (count > 0) {
      matched += 1;
      available.set(key, count - 1);
    }
  }
  return matched;
}

function detectParityCase(
  rdkit: RDKitModule,
  parityCase: ParityCase,
): readonly InteractionMetadata[] {
  const ligand = componentFromRdkitInput(
    rdkit,
    parityCase.ligand.molblock,
    parityCase.ligand.atoms,
    { removeHs: false },
  );
  const protein = componentFromRdkitInput(
    rdkit,
    parityCase.protein.molblock,
    parityCase.protein.atoms,
    { removeHs: false },
  );
  try {
    const factory = RULES[parityCase.rule];
    if (factory === undefined) {
      throw new Error(`Missing browser rule ${parityCase.rule}`);
    }
    return factory().detect(ligand, protein);
  } finally {
    ligand.dispose();
    protein.dispose();
  }
}

describe(`browser parity with ${fixture.generatedBy}`, () => {
  let rdkit: RDKitModule;

  beforeAll(async () => {
    rdkit = await loadRdkit();
  });

  it("only advertises rules covered by the Python oracle", () => {
    expect(fixture.rules).toEqual(Object.keys(RULES).sort());
    expect(fixture.pythonRdkitVersion).toMatch(/^2025\./);
    expect(rdkit.version()).toMatch(/^2025\./);
  });

  it("keeps positive and negative oracle coverage for every implemented rule", () => {
    expect(fixture.cases.length).toBeGreaterThanOrEqual(35);
    for (const rule of fixture.rules) {
      const cases = fixture.cases.filter((parityCase) => parityCase.rule === rule);
      expect(
        cases.some((parityCase) => parityCase.expected.length > 0),
        `${rule} is missing a positive oracle case`,
      ).toBe(true);
      expect(
        cases.some((parityCase) => parityCase.expected.length === 0),
        `${rule} is missing a negative oracle case`,
      ).toBe(true);
    }
  });

  it("measures 100% event precision and recall on the committed oracle corpus", () => {
    let expectedEvents = 0;
    let actualEvents = 0;
    let truePositives = 0;
    for (const parityCase of fixture.cases) {
      const actual = detectParityCase(rdkit, parityCase);
      expectedEvents += parityCase.expected.length;
      actualEvents += actual.length;
      truePositives += multisetIntersectionSize(actual, parityCase.expected);
    }
    expect(expectedEvents).toBeGreaterThanOrEqual(20);
    expect(actualEvents).toBeGreaterThanOrEqual(20);
    expect(truePositives / actualEvents, "event precision").toBe(1);
    expect(truePositives / expectedEvents, "event recall").toBe(1);
  });

  for (const parityCase of fixture.cases) {
    it(`matches atom participation and geometry: ${parityCase.name}`, () => {
      const actual = detectParityCase(rdkit, parityCase);
      expect(actual).toHaveLength(parityCase.expected.length);
      actual.forEach((interaction, index) => {
        const expected = parityCase.expected[index];
        expect(expected).toBeDefined();
        expect(interaction.interaction).toBe(expected?.interaction);
        expect(interaction.indices).toEqual(expected?.indices);
        expect(interaction.parentIndices).toEqual(expected?.parentIndices);
        expect(interaction.distance).toBeCloseTo(expected?.distance ?? NaN, 6);
        expect(Object.keys(interaction.geometry ?? {}).sort()).toEqual(
          Object.keys(expected?.geometry ?? {}).sort(),
        );
        for (const [key, value] of Object.entries(interaction.geometry ?? {})) {
          const expectedValue = expected?.geometry?.[key];
          if (Array.isArray(value)) {
            expect(Array.isArray(expectedValue)).toBe(true);
            expect(value).toHaveLength(
              Array.isArray(expectedValue) ? expectedValue.length : 0,
            );
            value.forEach((item, itemIndex) => {
              expect(item).toBeCloseTo(
                Array.isArray(expectedValue)
                  ? (expectedValue[itemIndex] ?? NaN)
                  : NaN,
                5,
              );
            });
          } else {
            expect(typeof expectedValue).toBe("number");
            expect(value).toBeCloseTo(
              typeof expectedValue === "number" ? expectedValue : NaN,
              5,
            );
          }
        }
      });
    });
  }

  it("keeps a broad Python-oracle SMARTS corpus", () => {
    expect(fixture.smartsCases.length).toBeGreaterThanOrEqual(150);
    expect(new Set(fixture.smartsCases.map((item) => item.patternKey))).toEqual(
      new Set(Object.keys(SMARTS_PATTERNS)),
    );
  });

  for (const smartsCase of fixture.smartsCases) {
    it(`matches Python SMARTS atom indices: ${smartsCase.name}`, () => {
      const patterns = SMARTS_PATTERNS[smartsCase.patternKey];
      expect(patterns, `missing browser pattern ${smartsCase.patternKey}`).toBeDefined();
      expect(smartsCase.pythonSmarts).toHaveLength(patterns?.length ?? 0);
      const component = componentFromRdkitInput(
        rdkit,
        smartsCase.component.molblock,
        smartsCase.component.atoms,
        { removeHs: false },
      );
      try {
        const actual = (patterns ?? []).flatMap((pattern) =>
          component.findMatches(pattern),
        );
        expect(actual).toEqual(smartsCase.expectedMatches);
      } finally {
        component.dispose();
      }
    });
  }
});
