import { describe, expect, it } from "vitest";

import { angleDegrees, centroid, distance, ringNormal } from "../src/index.js";

describe("geometry", () => {
  it("measures Euclidean distance in angstroms", () => {
    expect(distance([0, 0, 0], [2, 3, 6])).toBe(7);
  });

  it("measures a three-point angle", () => {
    expect(angleDegrees([1, 0, 0], [0, 0, 0], [0, 1, 0])).toBeCloseTo(90);
  });

  it("calculates centroids", () => {
    expect(centroid([[0, 0, 0], [3, 6, 9]])).toEqual([1.5, 3, 4.5]);
  });

  it("constructs a normalized ring-plane normal", () => {
    const normal = ringNormal([
      [1, 0, 0],
      [0, 1, 0],
      [-1, 0, 0],
      [0, -1, 0],
    ]);
    expect(Math.abs(normal[0])).toBeCloseTo(0);
    expect(Math.abs(normal[1])).toBeCloseTo(0);
    expect(Math.abs(normal[2])).toBeCloseTo(1);
  });

  it("rejects degenerate ring geometry", () => {
    expect(() => ringNormal([[0, 0, 0], [1, 0, 0], [2, 0, 0]])).toThrow(
      "zero-length vector",
    );
  });
});
