import { readFileSync } from "node:fs";

import type { RDKitLoader, RDKitModule } from "@rdkit/rdkit";
import { beforeAll, describe, expect, it } from "vitest";

import {
  componentFromRdkitInput,
  Hydrophobic,
  type AtomRecord,
  type InteractionMetadata,
} from "../src/index.js";

interface SerializedComponent {
  readonly molblock: string;
  readonly atoms: readonly AtomRecord[];
}

interface ParityCase {
  readonly name: string;
  readonly ligand: SerializedComponent;
  readonly protein: SerializedComponent;
  readonly expected: readonly InteractionMetadata[];
}

interface ParityFixture {
  readonly generatedBy: string;
  readonly pythonRdkitVersion: string;
  readonly rule: string;
  readonly cases: readonly ParityCase[];
}

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
    new URL("./reference/prolif-2.2.1-hydrophobic.json", import.meta.url),
    "utf8",
  ),
) as ParityFixture;

describe(`browser parity with ${fixture.generatedBy}`, () => {
  let rdkit: RDKitModule;

  beforeAll(async () => {
    rdkit = await loadRdkit();
  });

  it("only claims parity for the implemented rule", () => {
    expect(fixture.rule).toBe("Hydrophobic");
    expect(fixture.pythonRdkitVersion).toMatch(/^2025\./);
    expect(rdkit.version()).toMatch(/^2025\./);
  });

  for (const parityCase of fixture.cases) {
    it(`matches atom participation and geometry: ${parityCase.name}`, () => {
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
        const actual = new Hydrophobic().detect(ligand, protein);
        expect(actual).toHaveLength(parityCase.expected.length);
        actual.forEach((interaction, index) => {
          const expected = parityCase.expected[index];
          expect(expected).toBeDefined();
          expect(interaction.interaction).toBe(expected?.interaction);
          expect(interaction.indices).toEqual(expected?.indices);
          expect(interaction.parentIndices).toEqual(expected?.parentIndices);
          expect(interaction.distance).toBeCloseTo(expected?.distance ?? NaN, 6);
        });
      } finally {
        ligand.dispose();
        protein.dispose();
      }
    });
  }
});
