/** linear interpolate */
export const lerp = (value, sourceMin, sourceMax, targetMin, targetMax) =>
  targetMin +
  ((value - sourceMin) / (sourceMax - sourceMin || 1)) *
    (targetMax - targetMin);
