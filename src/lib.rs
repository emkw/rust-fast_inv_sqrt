//! Famous "Fast inverse square root" algorithm implementation
//! along with obligatory original comments.

#![cfg_attr(all(feature = "nightly",test), feature(asm,test,core_intrinsics))]

#[cfg(all(feature = "nightly",test))]
extern crate test;

#[cfg(test)]
#[macro_use] extern crate more_asserts;

use std::mem;
use std::f32;
use std::f64;

/// A trait that allows calculation of inverse square root into a Float(32).
pub trait InvSqrt32 {
	/// Calculates the inverse square root of a number
	fn inv_sqrt32(self) -> f32;
}

/// A trait that allows calculation of inverse square root into a Float(64).
pub trait InvSqrt64 {
	/// Calculates the inverse square root of a number.
	fn inv_sqrt64(self) -> f64;
}

/// Implementation of InvSqrt32 using the "Fast inverse square root" method.
/// This implementation somewhat provides look'n'feel of the original one.
/// If you want to uncomment 2nd iteration you actually have to make y mutable.
#[allow(non_upper_case_globals)]
impl InvSqrt32 for f32 {
	fn inv_sqrt32(self: f32) -> f32 {
		if cfg!(not(feature = "omit-checking")) {
			if self.signum() != 1.0 {
				return f32::NAN;
			} else if self == f32::INFINITY {
				return 0.0;
			} else if self < f32::MIN_POSITIVE {
				return f32::INFINITY;
			}
		}

		// Magic number based on Chris Lomont work:
		// const MAGIC_U32: u32 = 0x5f375a86;
		// The Original Magic Number:
		// const MAGIC_32: u32 = 0x5f3759df;
		const threehalfs: f32 = 1.5f32;
		let x2: f32 = self * 0.5f32;
		let mut i: u32 = unsafe { mem::transmute(self) }; // evil floating point bit level hacking
		i = 0x5f375a86 - (i >> 1);                        // what the fuck?
		let y: f32 = unsafe { mem::transmute(i) };
		let y  = y * ( threehalfs - ( x2 * y * y ) );     // 1st iteration
//		y  = y * ( threehalfs - ( x2 * y * y ) );       // 2nd iteration, this can be removed

		return y;
	}
}

/// Implementation of InvSqrt64 using the "Fast inverse square root" method.
/// This implementation provides somewhat more Rusty look'n'feel.
impl InvSqrt64 for f64 {
	fn inv_sqrt64(self: f64) -> f64 {
		if cfg!(not(feature = "omit-checking")) {
			if self.signum() != 1.0 {
				return f64::NAN;
			} else if self == f64::INFINITY {
				return 0.0;
			} else if self < f64::MIN_POSITIVE {
				return f64::INFINITY;
			}
		}

		// Magic number based on Chris Lomont work:
		const MAGIC_U64: u64 = 0x5fe6ec85e7de30da;
		const THREEHALFS: f64 = 1.5;
		let x2 = self * 0.5;
		let i = MAGIC_U64 - ( unsafe { mem::transmute::<_, u64>(self) } >> 1);
		let y: f64 = unsafe { mem::transmute(i) };

		y * ( THREEHALFS - ( x2 * y * y ) )
	}
}

impl InvSqrt32 for f64 {
	fn inv_sqrt32(self) -> f32 {
		(self as f32).inv_sqrt32()
	}
}

impl InvSqrt32 for i8 {
	fn inv_sqrt32(self) -> f32 {
		(self as f32).inv_sqrt32()
	}
}

impl InvSqrt32 for u8 {
	fn inv_sqrt32(self) -> f32 {
		(self as f32).inv_sqrt32()
	}
}

impl InvSqrt32 for i16 {
	fn inv_sqrt32(self) -> f32 {
		(self as f32).inv_sqrt32()
	}
}

impl InvSqrt32 for u16 {
	fn inv_sqrt32(self) -> f32 {
		(self as f32).inv_sqrt32()
	}
}

impl InvSqrt32 for i32 {
	fn inv_sqrt32(self) -> f32 {
		(self as f32).inv_sqrt32()
	}
}

impl InvSqrt32 for u32 {
	fn inv_sqrt32(self) -> f32 {
		(self as f32).inv_sqrt32()
	}
}

impl InvSqrt32 for i64 {
	fn inv_sqrt32(self) -> f32 {
		(self as f32).inv_sqrt32()
	}
}

impl InvSqrt32 for u64 {
	fn inv_sqrt32(self) -> f32 {
		(self as f32).inv_sqrt32()
	}
}

impl InvSqrt32 for isize {
	fn inv_sqrt32(self) -> f32 {
		(self as f32).inv_sqrt32()
	}
}

impl InvSqrt32 for usize {
	fn inv_sqrt32(self) -> f32 {
		(self as f32).inv_sqrt32()
	}
}

impl InvSqrt64 for f32 {
	fn inv_sqrt64(self) -> f64 {
		(self as f64).inv_sqrt64()
	}
}

impl InvSqrt64 for i8 {
	fn inv_sqrt64(self) -> f64 {
		(self as f64).inv_sqrt64()
	}
}

impl InvSqrt64 for u8 {
	fn inv_sqrt64(self) -> f64 {
		(self as f64).inv_sqrt64()
	}
}

impl InvSqrt64 for i16 {
	fn inv_sqrt64(self) -> f64 {
		(self as f64).inv_sqrt64()
	}
}

impl InvSqrt64 for u16 {
	fn inv_sqrt64(self) -> f64 {
		(self as f64).inv_sqrt64()
	}
}

impl InvSqrt64 for i32 {
	fn inv_sqrt64(self) -> f64 {
		(self as f64).inv_sqrt64()
	}
}

impl InvSqrt64 for u32 {
	fn inv_sqrt64(self) -> f64 {
		(self as f64).inv_sqrt64()
	}
}

impl InvSqrt64 for i64 {
	fn inv_sqrt64(self) -> f64 {
		(self as f64).inv_sqrt64()
	}
}

impl InvSqrt64 for u64 {
	fn inv_sqrt64(self) -> f64 {
		(self as f64).inv_sqrt64()
	}
}

impl InvSqrt64 for isize {
	fn inv_sqrt64(self) -> f64 {
		(self as f64).inv_sqrt64()
	}
}

impl InvSqrt64 for usize {
	fn inv_sqrt64(self) -> f64 {
		(self as f64).inv_sqrt64()
	}
}

#[cfg(test)]
mod test32 {
	use std::f32;
	use super::InvSqrt32;

	#[cfg(feature = "nightly")]
	use test::{black_box,Bencher};

	#[cfg(all(feature = "nightly", any(target_arch = "x86", target_arch = "x86_64")))]
	mod asm_x86 {
		use test::{black_box,Bencher};

		#[inline(always)]
		fn rsqrtss_mem(f: f32) -> f32 {
			let ret;
			unsafe {
				asm!("rsqrtss $1, $0" : "=x"(ret) : "m"(f) );
			}

			ret
		}

		#[inline(always)]
		fn rsqrtss_reg(f: f32) -> f32 {
			let ret;
			unsafe {
				asm!("rsqrtss $1, $0" : "=x"(ret) : "x"(f) );
			}

			ret
		}

		#[test]
		fn test_rsqrtss_mem() {
			let value: f32 = 11.1111;
			let result = rsqrtss_mem(value);
			assert_le!((0.3 - result).abs(), 0.0005);
		}

		#[test]
		fn test_rsqrtss_reg() {
			let value: f32 = 11.1111;
			let result = rsqrtss_reg(value);
			assert_le!((0.3 - result).abs(), 0.0005);
		}

		#[bench]
		fn bench_plain_rsqrtss_mem(b: &mut Bencher) {
			b.iter(|| {
				let f = black_box(1.2345f32);
				rsqrtss_mem(f)
			});
		}

		#[bench]
		fn bench_plain_rsqrtss_reg(b: &mut Bencher) {
			b.iter(|| {
				let f = black_box(1.2345f32);
				rsqrtss_reg(f)
			});
		}

		#[bench]
		fn bench_real_rsqrtss_reg_real(b: &mut Bencher) {
			b.iter(|| {
				let f = black_box(1.2345f32);
				let v = black_box(2.3456f32);
				v * rsqrtss_reg(f)
			});
		}
	}

	#[cfg(feature = "nightly")]
	mod llvm {
		use test::{black_box,Bencher};

		#[inline(always)]
		fn llvm_sqrt32(f: f32) -> f32 {
			unsafe { ::std::intrinsics::sqrtf32(f) }
		}

		#[inline(always)]
		fn llvm_inv_sqrt32(f: f32) -> f32 {
			1.0 / llvm_sqrt32(f)
		}

		#[test]
		fn test_llvm() {
			let value: f32 = 11.1111;
			let result = llvm_inv_sqrt32(value);
			assert_le!((0.3 - result).abs(), 0.0005);
		}

		#[bench]
		fn bench_plain_llvm(b: &mut Bencher) {
			b.iter(|| {
				let f = black_box(1.2345f32);
				llvm_inv_sqrt32(f)
			});
		}

		#[bench]
		fn bench_real_llvm(b: &mut Bencher) {
			b.iter(|| {
				let f = black_box(1.2345f32);
				let v = black_box(2.3456f32);
				v / llvm_sqrt32(f)
			});
		}
	}

	#[inline(always)]
	fn ref_inv_sqrt32(f: f32) -> f32 {
		1.0f32 / f.sqrt()
	}

	fn relative_difference(lhs: f32, rhs: f32) -> f32 {
		return 2.0 * (rhs - lhs).abs() / (rhs.abs() + lhs.abs());
	}

	#[test]
	fn test_f32() {
		let value: f32 = 11.1111;
		assert_le!((0.3 - value.inv_sqrt32()).abs(), 0.0005);

		let mut value: f32 = 0.00000001f32;
		while value < 100.4f32 {
			let result = value.inv_sqrt32();
			let ref_result = ref_inv_sqrt32(value);
			let relative_diff = relative_difference(result, ref_result);
			//println!("inv_sqrt32() -> {}: {} / {} [{}]", value, result, ref_result, relative_diff);
			assert_le!(relative_diff, 0.002);
			value += 0.01f32;
		}
	}

	#[test]
	fn test_zero() {
		let zero = 0.0f32;
		assert_eq!(zero.inv_sqrt32(), f32::INFINITY);
	}

	#[test]
	fn test_negative() {
		let negative = -1.0f32;
		assert!(negative.inv_sqrt32().is_nan());
	}

	#[test]
	fn test_negative_zero() {
		let negative_zero = -0.0f32;
		assert!(negative_zero.inv_sqrt32().is_nan());
	}

	#[test]
	fn test_i8() {
		let value: i8 = 55;
		assert_le!((0.1348399725 - value.inv_sqrt32()).abs(), 0.00001);
	}

	#[cfg(feature = "nightly")]
	#[bench]
	fn bench_plain_ref(b: &mut Bencher) {
		b.iter(|| {
			let f = black_box(1.2345f32);
			ref_inv_sqrt32(f)
		});
	}

	#[cfg(feature = "nightly")]
	#[bench]
	fn bench_recip_sqrt(b: &mut Bencher) {
		b.iter(|| {
			let f = black_box(1.2345f32);
			f.recip().sqrt()
		});
	}

	#[cfg(feature = "nightly")]
	#[bench]
	fn bench_sqrt_recip(b: &mut Bencher) {
		b.iter(|| {
			let f = black_box(1.2345f32);
			f.sqrt().recip()
		});
	}

	#[cfg(feature = "nightly")]
	#[bench]
	fn bench_plain_impl(b: &mut Bencher) {
		b.iter(|| {
			let f = black_box(1.2345f32);
			f.inv_sqrt32()
		});
	}

	#[cfg(feature = "nightly")]
	#[bench]
	fn bench_real_ref(b: &mut Bencher) {
		b.iter(|| {
			let v = black_box(2.3456f32);
			let f = black_box(1.2345f32);
			v / f.sqrt()
		});
	}

	#[cfg(feature = "nightly")]
	#[bench]
	fn bench_real_impl(b: &mut Bencher) {
		b.iter(|| {
			let v = black_box(2.3456f32);
			let f = black_box(1.2345f32);
			v * f.inv_sqrt32()
		});
	}
}

#[cfg(test)]
mod test64 {
	use std::f64;
	use super::InvSqrt64;

	#[cfg(feature = "nightly")]
	use test::{black_box,Bencher};

	#[inline(always)]
	fn ref_inv_sqrt64(f: f64) -> f64 {
		1.0f64/f.sqrt()
	}

	fn relative_difference(lhs: f64, rhs: f64) -> f64 {
		return 2.0 * (rhs - lhs).abs() / (rhs.abs() + lhs.abs());
	}

	#[test]
	fn test_f64() {
		let value: f64 = 11.1111;
		assert_le!((0.3 - value.inv_sqrt64()).abs(), 0.0005);

		let mut value: f64 = 0.00000001f64;
		while value < 100.4f64 {
			let result = value.inv_sqrt64();
			let ref_result = ref_inv_sqrt64(value);
			let relative_diff = relative_difference(result, ref_result);
			//println!("inv_sqrt64(): {} -> {} / {} [{}]", value, result, ref_result, relative_diff);
			assert_le!(relative_diff, 0.002);
			value += 0.01f64;
		}
	}

	#[test]
	fn test_zero() {
		let zero = 0.0f64;
		assert_eq!(zero.inv_sqrt64(), f64::INFINITY);
	}

	#[test]
	fn test_negative() {
		let negative = -1.0f64;
		assert!(negative.inv_sqrt64().is_nan());
	}

	#[test]
	fn test_negative_zero() {
		let negative_zero = -0.0f64;
		assert!(negative_zero.inv_sqrt64().is_nan());
	}

	#[test]
	fn test_i8() {
		let value: i8 = 55;
		assert_le!((0.1348399725 - value.inv_sqrt64()).abs(), 0.00001);
	}

	#[cfg(feature = "nightly")]
	#[bench]
	fn bench_plain_ref(b: &mut Bencher) {
		b.iter(|| {
			let f = black_box(1.2345f64);
			ref_inv_sqrt64(f)
		});
	}

	#[cfg(feature = "nightly")]
	#[bench]
	fn bench_recip_sqrt(b: &mut Bencher) {
		b.iter(|| {
			let f = black_box(1.2345f64);
			f.recip().sqrt()
		});
	}

	#[cfg(feature = "nightly")]
	#[bench]
	fn bench_sqrt_recip(b: &mut Bencher) {
		b.iter(|| {
			let f = black_box(1.2345f64);
			f.sqrt().recip()
		});
	}

	#[cfg(feature = "nightly")]
	#[bench]
	fn bench_plain_impl(b: &mut Bencher) {
		b.iter(|| {
			let f = black_box(1.2345f64);
			f.inv_sqrt64()
		});
	}

	#[cfg(feature = "nightly")]
	#[bench]
	fn bench_real_ref(b: &mut Bencher) {
		b.iter(|| {
			let v = black_box(2.3456f64);
			let f = black_box(1.2345f64);
			v / f.sqrt()
		});
	}

	#[cfg(feature = "nightly")]
	#[bench]
	fn bench_real_impl(b: &mut Bencher) {
		b.iter(|| {
			let v = black_box(2.3456f64);
			let f = black_box(1.2345f64);
			v * f.inv_sqrt64()
		});
	}
}
