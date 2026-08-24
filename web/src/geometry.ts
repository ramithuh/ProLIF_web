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

export function add(a: Vec3, b: Vec3): Vec3 {
  return [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
}

export function scale(vector: Vec3, factor: number): Vec3 {
  return [vector[0] * factor, vector[1] * factor, vector[2] * factor];
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

export function vectorAngleDegrees(first: Vec3, second: Vec3): number {
  const a = normalize(first);
  const b = normalize(second);
  const cosine = Math.min(1, Math.max(-1, dot(a, b)));
  return (Math.acos(cosine) * 180) / Math.PI;
}

/** Fold a directed normal-vector angle onto the unoriented 0–90° ring range. */
export function foldRingAngle(angle: number): number {
  if (angle >= 180) return angle % 90;
  if (angle > 90) return 90 - (angle % 90);
  return angle;
}

export function betweenInclusive(
  value: number,
  limits: readonly [minimum: number, maximum: number],
): boolean {
  return value >= limits[0] && value <= limits[1];
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

function determinant3(matrix: readonly [Vec3, Vec3, Vec3]): number {
  const [a, b, c] = matrix;
  return (
    a[0] * (b[1] * c[2] - b[2] * c[1]) -
    a[1] * (b[0] * c[2] - b[2] * c[0]) +
    a[2] * (b[0] * c[1] - b[1] * c[0])
  );
}

/** Solve a 3×3 linear system using Cramer's rule. */
export function solve3(
  matrix: readonly [Vec3, Vec3, Vec3],
  values: Vec3,
): Vec3 | null {
  const determinant = determinant3(matrix);
  if (Math.abs(determinant) <= EPSILON) return null;
  const [row0, row1, row2] = matrix;
  const column0: Vec3 = [values[0], values[1], values[2]];
  const column1: Vec3 = [row0[1], row1[1], row2[1]];
  const column2: Vec3 = [row0[2], row1[2], row2[2]];
  const original0: Vec3 = [row0[0], row1[0], row2[0]];

  const x = determinant3([
    [column0[0], column1[0], column2[0]],
    [column0[1], column1[1], column2[1]],
    [column0[2], column1[2], column2[2]],
  ]);
  const y = determinant3([
    [original0[0], column0[0], column2[0]],
    [original0[1], column0[1], column2[1]],
    [original0[2], column0[2], column2[2]],
  ]);
  const z = determinant3([
    [original0[0], column1[0], column0[0]],
    [original0[1], column1[1], column0[1]],
    [original0[2], column1[2], column0[2]],
  ]);
  return [x / determinant, y / determinant, z / determinant];
}

/** Match ProLIF's line-of-intersection projection used by edge-to-face rings. */
export function ringPlaneIntersection(
  planeNormal: Vec3,
  planeCentroid: Vec3,
  tiltedNormal: Vec3,
  tiltedCentroid: Vec3,
): Vec3 | null {
  const direction = cross(planeNormal, tiltedNormal);
  const point = solve3(
    [planeNormal, tiltedNormal, direction],
    [dot(planeNormal, planeCentroid), dot(tiltedNormal, tiltedCentroid), 0],
  );
  if (point === null) return null;
  const unitDirection = normalize(direction);
  const scalarProjection = dot(unitDirection, subtract(planeCentroid, point));
  return add(point, scale(unitDirection, scalarProjection));
}
