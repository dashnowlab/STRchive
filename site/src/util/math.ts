/** linear interpolate */
export const lerp = (
  value: number,
  sourceMin: number,
  sourceMax: number,
  targetMin: number,
  targetMax: number,
) =>
  targetMin +
  ((value - sourceMin) / (sourceMax - sourceMin || 1)) *
    (targetMax - targetMin);
