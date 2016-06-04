# fast_inv_sqrt

This was made purely for fun and testing crates.io publishing, but may actually be usable.

InvSqrt32 trait provide inv_sqrt32() function for primitive numeric types.
InvSqrt64 provides inv_sqrt64().

## How do I use it?

Cargo.toml:

```toml
[dependencies]
fast_inv_sqrt = "1.0.0"
```

Code:

```rust
extern crate fast_inv_sqrt;

use fast_inv_sqrt::InvSqrt32;
use fast_inv_sqrt::InvSqrt64;

fn main() {
	let f: f32 = 1.234;
	println!("{}", f.inv_sqrt32());

	let i: i8 = 55;
	println!("{}", i.inv_sqrt64());
}
```

## Benchmarks

Benchmarks require nightly compiler.

"_ref_" benchmarks use f1 / f2.sqrt()
"_impl_" benchmarks use f1 * f2.inv_sqrtXX()

$ cargo bench --features 'nightly'
```
test test32::bench_plain_impl ... bench:           8 ns/iter (+/- 0)
test test32::bench_plain_ref  ... bench:          13 ns/iter (+/- 0)

test test32::bench_real_impl  ... bench:           8 ns/iter (+/- 0)
test test32::bench_real_ref   ... bench:          13 ns/iter (+/- 0)

test test64::bench_plain_impl ... bench:          10 ns/iter (+/- 0)
test test64::bench_plain_ref  ... bench:          20 ns/iter (+/- 0)

test test64::bench_real_impl  ... bench:          10 ns/iter (+/- 0)
test test64::bench_real_ref   ... bench:          20 ns/iter (+/- 0)
```
Feature 'omit-checking' disables checks of if value passed is_sign_positive() and 
is_normal(), and will produce invalid results for denormalized and negative values,
but this is fast:
```
test test32::bench_plain_impl ... bench:           2 ns/iter (+/- 0)
test test32::bench_plain_ref  ... bench:          13 ns/iter (+/- 0)

test test32::bench_real_impl  ... bench:           3 ns/iter (+/- 0)
test test32::bench_real_ref   ... bench:          13 ns/iter (+/- 0)

test test64::bench_plain_impl ... bench:           2 ns/iter (+/- 0)
test test64::bench_plain_ref  ... bench:          20 ns/iter (+/- 0)

test test64::bench_real_impl  ... bench:           3 ns/iter (+/- 0)
test test64::bench_real_ref   ... bench:          20 ns/iter (+/- 0)
```
