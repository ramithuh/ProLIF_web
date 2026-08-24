import type { Vec3 } from "./types.js";

const EPSILON = 1e-12;

export function subtract(a: Vec3, b: Vec3): Vec3 {
  return [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
}

export function dot(a: Vec3, b: Vec3): number {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

export function cross(a: Vec3, b: Vec3): Vec3 {
  return [
    a[1] * b[2] - a[2] * b[1],
    a[2] * b[0] - a[0] * b[2],
    a[0] * b[1] - a[1] * b[0],
  ];
}

export function norm(vector: Vec3): number {
  return Math.sqrt(dot(vector, vector));
}

export function normalize(vector: Vec3): Vec3 {
  const length = norm(vector);
  if (length <= EPSILON) {
    throw new RangeError("Cannot normalize a zero-length vector");
  }
  return [vector[0] / length, vector[1] / length, vector[2] / length];
}

export function distance(a: Vec3, b: Vec3): number {
  return norm(subtract(a, b));
}

export function angleDegrees(a: Vec3, vertex: Vec3, c: Vec3): number {
  const first = normalize(subtract(a, vertex));
  const second = normalize(subtract(c, vertex));
  const cosine = Math.min(1, Math.max(-1, dot(first, second)));
  return (Math.acos(cosine) * 180) / Math.PI;
}

export function centroid(points: readonly Vec3[]): Vec3 {
  if (points.length === 0) {
    throw new RangeError("Cannot calculate the centroid of an empty point set");
  }
  const sum = points.reduce<Vec3>(
    (total, point) => [
      total[0] + point[0],
      total[1] + point[1],
      total[2] + point[2],
    ],
    [0, 0, 0],
  );
  return [sum[0] / points.length, sum[1] / points.length, sum[2] / points.length];
}

/**
 * Unit normal for an ordered ring. This follows ProLIF's centroid-to-first-two
 * atoms construction and deliberately rejects degenerate rings.
 */
export function ringNormal(points: readonly Vec3[]): Vec3 {
  if (points.length < 3) {
    throw new RangeError("At least three points are required for a ring normal");
  }
  const center = centroid(points);
  const first = points[0];
  const second = points[1];
  if (first === undefined || second === undefined) {
    throw new RangeError("Ring coordinates are incomplete");
  }
  return normalize(cross(subtract(first, center), subtract(second, center)));
}
