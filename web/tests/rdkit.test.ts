import type { RDKitLoader, RDKitModule } from "@rdkit/rdkit";
import { afterAll, beforeAll, describe, expect, it } from "vitest";

import {
  componentFromRdkitInput,
  PROLIF_2_2_HYDROPHOBIC_SMARTS,
} from "../src/index.js";
import { atom } from "./fixtures.js";

async function loadRdkit(): Promise<RDKitModule> {
  const namespace: unknown = await import("@rdkit/rdkit");
  const candidate = namespace as {
    readonly default?: RDKitLoader;
  };
  if (typeof candidate.default !== "function") {
    throw new TypeError("@rdkit/rdkit did not expose its module loader");
  }
  return candidate.default();
}

describe("RDKit.js SMARTS adapter", () => {
  let rdkit: RDKitModule;

  beforeAll(async () => {
    rdkit = await loadRdkit();
  });

  afterAll(() => {
    // The module itself has no public disposal method; owned molecules are
    // released by each component below.
  });

  it.each([
    ["C", 0],
    ["c1cscc1", 5],
    ["c1cocc1", 2],
    ["BrI", 2],
    ["C=O", 0],
    ["Nc1ccccc1", 5],
  ])(
    "matches ProLIF's hydrophobic SMARTS count for %s",
    (smiles, expectedCount) => {
      const molecule = rdkit.get_mol(smiles);
      if (molecule === null) {
        throw new Error(`Could not parse test molecule ${smiles}`);
      }
      const moleculeJson = JSON.parse(molecule.get_json()) as {
        molecules: readonly { atoms: readonly unknown[] }[];
      };
      const atomCount = moleculeJson.molecules[0]?.atoms.length;
      molecule.delete();
      if (atomCount === undefined) {
        throw new Error(`Could not read atom count for ${smiles}`);
      }

      const component = componentFromRdkitInput(
        rdkit,
        smiles,
        Array.from({ length: atomCount }, (_, index) =>
          atom(index, index, [index, 0, 0]),
        ),
      );
      try {
        expect(
          component.findMatches(PROLIF_2_2_HYDROPHOBIC_SMARTS),
        ).toHaveLength(expectedCount);
      } finally {
        component.dispose();
      }
    },
  );

  it("rejects invalid SMARTS without leaking the molecule", () => {
    const component = componentFromRdkitInput(rdkit, "CC", [
      atom(0, 0, [0, 0, 0]),
      atom(1, 1, [1, 0, 0]),
    ]);
    try {
      expect(() => component.findMatches("[not-valid")).toThrow(
        "could not parse SMARTS",
      );
    } finally {
      component.dispose();
    }
  });

  it("rejects coordinate arrays that do not follow RDKit atom order", () => {
    expect(() =>
      componentFromRdkitInput(rdkit, "C", [atom(7, 7, [0, 0, 0])]),
    ).toThrow("must follow RDKit ordering");
  });
});
