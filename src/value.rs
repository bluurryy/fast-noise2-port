use crate::floor;

use super::{interp_hermite, lerp, primes, value_coord2, value_coord3, value_coord4};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Value;

impl Value {
    #[inline]
    pub fn sample2(self, [x, y]: [f32; 2], seed: i32) -> f32 {
        let xs = floor(x);
        let ys = floor(y);

        let x0 = (xs as i32).wrapping_mul(primes::X);
        let y0 = (ys as i32).wrapping_mul(primes::Y);

        let x1 = x0.wrapping_add(primes::X);
        let y1 = y0.wrapping_add(primes::Y);

        let xs = interp_hermite(x - xs);
        let ys = interp_hermite(y - ys);

        lerp(
            lerp(value_coord2(seed, x0, y0), value_coord2(seed, x1, y0), xs),
            lerp(value_coord2(seed, x0, y1), value_coord2(seed, x1, y1), xs),
            ys,
        )
    }

    #[inline]
    pub fn sample3(self, [x, y, z]: [f32; 3], seed: i32) -> f32 {
        let xs = floor(x);
        let ys = floor(y);
        let zs = floor(z);

        let x0 = (xs as i32).wrapping_mul(primes::X);
        let y0 = (ys as i32).wrapping_mul(primes::Y);
        let z0 = (zs as i32).wrapping_mul(primes::Z);

        let x1 = x0.wrapping_add(primes::X);
        let y1 = y0.wrapping_add(primes::Y);
        let z1 = z0.wrapping_add(primes::Z);

        let xs = interp_hermite(x - xs);
        let ys = interp_hermite(y - ys);
        let zs = interp_hermite(z - zs);

        lerp(
            lerp(
                lerp(
                    value_coord3(seed, x0, y0, z0),
                    value_coord3(seed, x1, y0, z0),
                    xs,
                ),
                lerp(
                    value_coord3(seed, x0, y1, z0),
                    value_coord3(seed, x1, y1, z0),
                    xs,
                ),
                ys,
            ),
            lerp(
                lerp(
                    value_coord3(seed, x0, y0, z1),
                    value_coord3(seed, x1, y0, z1),
                    xs,
                ),
                lerp(
                    value_coord3(seed, x0, y1, z1),
                    value_coord3(seed, x1, y1, z1),
                    xs,
                ),
                ys,
            ),
            zs,
        )
    }

    #[inline]
    pub fn sample4(self, [x, y, z, w]: [f32; 4], seed: i32) -> f32 {
        let xs = floor(x);
        let ys = floor(y);
        let zs = floor(z);
        let ws = floor(w);

        let x0 = (xs as i32).wrapping_mul(primes::X);
        let y0 = (ys as i32).wrapping_mul(primes::Y);
        let z0 = (zs as i32).wrapping_mul(primes::Z);
        let w0 = (ws as i32).wrapping_mul(primes::W);

        let x1 = x0.wrapping_add(primes::X);
        let y1 = y0.wrapping_add(primes::Y);
        let z1 = z0.wrapping_add(primes::Z);
        let w1 = w0.wrapping_add(primes::W);

        let xs = interp_hermite(x - xs);
        let ys = interp_hermite(y - ys);
        let zs = interp_hermite(z - zs);
        let ws = interp_hermite(w - ws);

        lerp(
            lerp(
                lerp(
                    lerp(
                        value_coord4(seed, x0, y0, z0, w0),
                        value_coord4(seed, x1, y0, z0, w0),
                        xs,
                    ),
                    lerp(
                        value_coord4(seed, x0, y1, z0, w0),
                        value_coord4(seed, x1, y1, z0, w0),
                        xs,
                    ),
                    ys,
                ),
                lerp(
                    lerp(
                        value_coord4(seed, x0, y0, z1, w0),
                        value_coord4(seed, x1, y0, z1, w0),
                        xs,
                    ),
                    lerp(
                        value_coord4(seed, x0, y1, z1, w0),
                        value_coord4(seed, x1, y1, z1, w0),
                        xs,
                    ),
                    ys,
                ),
                zs,
            ),
            lerp(
                lerp(
                    lerp(
                        value_coord4(seed, x0, y0, z0, w1),
                        value_coord4(seed, x1, y0, z0, w1),
                        xs,
                    ),
                    lerp(
                        value_coord4(seed, x0, y1, z0, w1),
                        value_coord4(seed, x1, y1, z0, w1),
                        xs,
                    ),
                    ys,
                ),
                lerp(
                    lerp(
                        value_coord4(seed, x0, y0, z1, w1),
                        value_coord4(seed, x1, y0, z1, w1),
                        xs,
                    ),
                    lerp(
                        value_coord4(seed, x0, y1, z1, w1),
                        value_coord4(seed, x1, y1, z1, w1),
                        xs,
                    ),
                    ys,
                ),
                zs,
            ),
            ws,
        )
    }
}
