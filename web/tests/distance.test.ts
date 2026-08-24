import { describe, expect, it } from "vitest";

import { DistanceRule, Hydrophobic } from "../src/index.js";
import { atom, FixtureComponent } from "./fixtures.js";

describe("DistanceRule", () => {
  it("matches ProLIF's inclusive first-match-atom cutoff behavior", () => {
    const ligand = new FixtureComponent(
      [atom(0, 10, [0, 0, 0]), atom(1, 11, [100, 0, 0])],
      { "lig-query": [[0, 1]] },
    );
    const protein = new FixtureComponent(
      [atom(0, 42, [4.5, 0, 0]), atom(1, 43, [100, 0, 0])],
      { "prot-query": [[0, 1]] },
    );
    const rule = new DistanceRule({
      name: "ReferenceDistance",
      ligandSmarts: "lig-query",
      proteinSmarts: "prot-query",
      distance: 4.5,
    });

    expect(rule.detect(ligand, protein)).toEqual([
      {
        interaction: "ReferenceDistance",
        indices: { ligand: [0, 1], protein: [0, 1] },
        parentIndices: { ligand: [10, 11], protein: [42, 43] },
        distance: 4.5,
      },
    ]);
  });

  it("returns every matching feature pair inside the cutoff", () => {
    const ligand = new FixtureComponent(
      [atom(0, 0, [0, 0, 0]), atom(1, 1, [10, 0, 0])],
      { q: [[0], [1]] },
    );
    const protein = new FixtureComponent(
      [atom(0, 2, [1, 0, 0]), atom(1, 3, [12, 0, 0])],
      { q: [[0], [1]] },
    );
    const rule = new DistanceRule({
      name: "AllPairs",
      ligandSmarts: "q",
      proteinSmarts: "q",
      distance: 3,
    });

    expect(rule.detect(ligand, protein).map((item) => item.distance)).toEqual([
      1, 2,
    ]);
  });

  it("rejects malformed matches instead of silently inventing atom indices", () => {
    const ligand = new FixtureComponent([atom(0, 0, [0, 0, 0])], { q: [[]] });
    const protein = new FixtureComponent([atom(0, 0, [0, 0, 0])], { q: [[0]] });
    const rule = new DistanceRule({
      name: "Malformed",
      ligandSmarts: "q",
      proteinSmarts: "q",
      distance: 1,
    });

    expect(() => rule.detect(ligand, protein)).toThrow("SMARTS match is empty");
  });
});

describe("Hydrophobic", () => {
  it("uses ProLIF 2.2.x's default 4.5 angstrom cutoff", () => {
    const rule = new Hydrophobic();
    expect(rule.cutoff).toBe(4.5);
    expect(rule.name).toBe("Hydrophobic");
    expect(rule.ligandSmarts).toBe(rule.proteinSmarts);
  });
});
