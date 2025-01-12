#![expect(dead_code)]

#[cfg(feature = "libm")]
mod libm_math {
    #[inline(always)]
    pub(crate) fn abs(f: f32) -> f32 {
        libm::fabsf(f)
    }

    #[inline(always)]
    pub(crate) fn sqrt(f: f32) -> f32 {
        libm::sqrtf(f)
    }

    #[inline(always)]
    pub(crate) fn trunc(f: f32) -> f32 {
        libm::truncf(f)
    }

    #[inline(always)]
    pub(crate) fn floor(f: f32) -> f32 {
        libm::floorf(f)
    }

    #[inline(always)]
    pub(crate) fn round(f: f32) -> f32 {
        libm::roundf(f)
    }

    #[inline(always)]
    pub(crate) fn sin(f: f32) -> f32 {
        libm::sinf(f)
    }

    #[inline(always)]
    pub(crate) fn cos(f: f32) -> f32 {
        libm::cosf(f)
    }

    #[inline(always)]
    pub(crate) fn min(a: f32, b: f32) -> f32 {
        libm::fminf(a, b)
    }

    #[inline(always)]
    pub(crate) fn max(a: f32, b: f32) -> f32 {
        libm::fmaxf(a, b)
    }

    #[inline(always)]
    pub(crate) fn mul_add(a: f32, b: f32, c: f32) -> f32 {
        a * b + c
    }
}

#[cfg(not(feature = "libm"))]
mod std_math {
    #[inline(always)]
    pub(crate) fn abs(f: f32) -> f32 {
        f32::abs(f)
    }

    #[inline(always)]
    pub(crate) fn sqrt(f: f32) -> f32 {
        f32::sqrt(f)
    }

    #[inline(always)]
    pub(crate) fn trunc(f: f32) -> f32 {
        f32::trunc(f)
    }

    #[inline(always)]
    pub(crate) fn floor(f: f32) -> f32 {
        f32::floor(f)
    }

    #[inline(always)]
    pub(crate) fn round(f: f32) -> f32 {
        f32::round(f)
    }

    #[inline(always)]
    pub(crate) fn sin(f: f32) -> f32 {
        f32::sin(f)
    }

    #[inline(always)]
    pub(crate) fn cos(f: f32) -> f32 {
        f32::cos(f)
    }

    #[inline(always)]
    pub(crate) fn min(a: f32, b: f32) -> f32 {
        f32::min(a, b)
    }

    #[inline(always)]
    pub(crate) fn max(a: f32, b: f32) -> f32 {
        f32::max(a, b)
    }

    #[inline(always)]
    pub(crate) fn mul_add(a: f32, b: f32, c: f32) -> f32 {
        f32::mul_add(a, b, c)
    }
}

#[cfg(feature = "libm")]
pub(crate) use libm_math::*;

#[cfg(not(feature = "libm"))]
pub(crate) use std_math::*;

pub trait Dot {
    type Output;
    fn dot(lhs: Self, rhs: Self) -> Self::Output;
}

#[inline(always)]
pub fn dot<T: Dot>(lhs: T, rhs: T) -> T::Output {
    Dot::dot(lhs, rhs)
}

#[inline(always)]
pub fn length_squared<T: Dot + Copy>(value: T) -> T::Output {
    dot(value, value)
}

#[allow(dead_code)]
pub fn length<T: Dot<Output = f32> + Copy>(value: T) -> T::Output {
    crate::math::sqrt(length_squared(value))
}

pub trait FloorToInt {
    type Output;
    fn floor_to_int(self) -> Self::Output;
}

impl FloorToInt for f32 {
    type Output = i32;

    #[inline(always)]
    fn floor_to_int(self) -> Self::Output {
        if self >= 0.0 {
            self as i32
        } else {
            (self as i32).wrapping_sub(1)
        }
    }
}

#[inline(always)]
pub fn floor_to_int<T: FloorToInt>(value: T) -> T::Output {
    FloorToInt::floor_to_int(value)
}

pub trait RoundToInt {
    type Output;
    fn round_to_int(self) -> Self::Output;
}

impl RoundToInt for f32 {
    type Output = i32;

    #[inline(always)]
    fn round_to_int(self) -> Self::Output {
        if self >= 0.0 {
            (self + 0.5) as i32
        } else {
            (self - 0.5) as i32
        }
    }
}

#[inline(always)]
pub fn round_to_int<T: RoundToInt>(value: T) -> T::Output {
    RoundToInt::round_to_int(value)
}

pub trait InterpQuintic {
    fn interp_quintic(self) -> Self;
}

impl InterpQuintic for f32 {
    #[inline(always)]
    fn interp_quintic(self) -> Self {
        self * self * self * (self * (self * 6.0 - 15.0) + 10.0)
    }
}

#[inline(always)]
pub fn interp_quintic<T: InterpQuintic>(value: T) -> T {
    InterpQuintic::interp_quintic(value)
}

pub trait InterpHermite {
    fn interp_hermite(self) -> Self;
}

impl InterpHermite for f32 {
    #[inline(always)]
    fn interp_hermite(self) -> Self {
        self * self * (3.0 - 2.0 * self)
    }
}

#[inline(always)]
pub fn interp_hermite<T: InterpHermite>(value: T) -> T {
    InterpHermite::interp_hermite(value)
}

pub trait Lerp<T = Self> {
    type Output;
    fn lerp(a: Self, b: Self, t: T) -> Self::Output;
}

impl Lerp for f32 {
    type Output = Self;

    #[inline(always)]
    fn lerp(a: Self, b: Self, t: Self) -> Self::Output {
        a + t * (b - a)
    }
}

#[inline(always)]
pub fn lerp<T, V: Lerp<T>>(a: V, b: V, t: T) -> V::Output {
    Lerp::lerp(a, b, t)
}
