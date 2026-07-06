# Splines

[![CI](https://github.com/pbosetti/Splines/actions/workflows/ci.yml/badge.svg)](https://github.com/pbosetti/Splines/actions/workflows/ci.yml)
![C++20](https://img.shields.io/badge/C%2B%2B-20-blue.svg)
![CMake](https://img.shields.io/badge/CMake-3.25%2B-informational.svg)

**A C++ library for univariate and bivariate spline interpolation** — linear,
Akima, Bessel, cubic, Hermite, PCHIP, and quintic curves, plus bilinear,
bicubic, and biquintic surfaces, all built from raw arrays, `nlohmann::json`,
YAML/JSON via `GenericContainer`, or a plain C interface for non-C++ hosts.

```cpp
#include "Splines/Splines.hh"

std::vector<double> x = { 0, 1, 2, 3, 4 };
std::vector<double> y = { 0, 1, 4, 9, 16 };

Splines::CubicSpline spline;
spline.build( x, y );

double value      = spline( 2.5 );   // interpolated value
double derivative = spline.D( 2.5 ); // first derivative

// Batch evaluation over an array of abscissas (eval / D / DD):
std::vector<double> xs = { 0.5, 1.5, 2.5, 3.5 }, ys( xs.size() );
spline.eval( xs, ys );               // ys[i] == spline(xs[i]), bit-for-bit
```

Batch results are identical to calling the scalar form per element; for
**sorted** abscissas (resampling, plotting) the batch path reuses the
current segment instead of re-searching, ~3-4x faster than the scalar loop.

---

## Contents

- [Splines](#splines)
  - [Contents](#contents)
  - [Curve and surface types](#curve-and-surface-types)
  - [Building](#building)
  - [Using it from your own CMake project](#using-it-from-your-own-cmake-project)
    - [`FetchContent` (recommended)](#fetchcontent-recommended)
    - [If other subprojects also fetch Splines](#if-other-subprojects-also-fetch-splines)
    - [Installing without CMake integration](#installing-without-cmake-integration)
  - [The C interface](#the-c-interface)
  - [Performance](#performance)
  - [Testing](#testing)
  - [License](#license)

---

## Curve and surface types

- Univariate (`Splines::Spline` subclasses): Linear, Akima, Bessel, Cubic,
  Hermite, PCHIP (non-oscillatory), Quintic
- Bivariate (`Splines::SplineSurf` subclasses): Bilinear, Bicubic, Biquintic
- `Splines::SplineSet` groups several 1D splines that share the same
  independent variable (e.g. a parametric path), building/evaluating them
  together and bridging to/from JSON and YAML via
  [`GenericContainer`](https://github.com/pbosetti/GenericContainer)

## Building

Requires **CMake ≥ 3.25** and a **C++20** compiler (Clang, GCC, or MSVC).
All dependencies are fetched automatically via `FetchContent`:
[nlohmann/json](https://github.com/nlohmann/json),
[GenericContainer](https://github.com/pbosetti/GenericContainer),
[UtilsLite](https://github.com/pbosetti/UtilsLite), and
[quarticRootsFlocke](https://github.com/ebertolazzi/quarticRootsFlocke).

> **Building with Clang on Linux**: `UtilsLite` and `quarticRootsFlocke`
> both build against `libc++` when the compiler is Clang, so Splines
> matches that choice to keep the final link ABI-consistent. macOS ships
> `libc++` with the system Clang; on Linux install it explicitly first
> (Debian/Ubuntu: `apt install libc++-dev libc++abi-dev`). Building with
> GCC instead sidesteps this entirely.

```sh
cmake -S . -B build                             # configure (Release by default)
cmake --build build                             # build
ctest --test-dir build                          # run the test suite
cmake --install build --prefix /path/to/prefix  # install
cmake --build build --target package            # cpack: .tar.gz / .zip
```

> **NOTE**: `build.sh` (macOS/Linux) and `build.ps1` (Windows) wrap the
> commands above into a single call, e.g. `./build.sh test Release`.

Useful configure-time options:

| Option | Default | Effect |
|---|---|---|
| `BUILD_SHARED_LIBS` | `OFF` | Build a shared (`.so`/`.dylib`/`.dll`) instead of a static library |
| `BUILD_TESTING` | `ON` | Build the `test00`..`test18` executables |
| `SPLINES_INSTALL` | `ON` (top-level) | Generate install/package rules |
| `SPLINES_BUILD_BENCHMARKS` | `ON` (top-level) | Build `bin/bench_eval` (micro-benchmark, run manually) |
| `SPLINES_STRICT_WARNINGS` | `OFF` | Compile the Splines target with `-Wconversion`/`-Wsign-conversion`/`-Wshadow`/`-Wdouble-promotion` (GCC/Clang) |

## Using it from your own CMake project

### `FetchContent` (recommended)

```cmake
include( FetchContent )
FetchContent_Declare(
  Splines
  GIT_REPOSITORY https://github.com/ebertolazzi/Splines.git
  GIT_TAG        develop   # pin to a commit or tag for reproducible builds
)
set( BUILD_TESTING OFF )   # don't build Splines' own test suite as a dependency
FetchContent_MakeAvailable( Splines )

add_executable( my_app main.cc )
target_link_libraries( my_app PRIVATE Splines::Splines )
```

`Splines::Splines` publicly links `GenericContainer::Yaml`, `UtilsLite`, and
`quarticRootsFlocke`, so their headers and link requirements reach `my_app`
transitively — no need to reference them directly unless you use them
yourself.

### If other subprojects also fetch Splines

`FetchContent` deduplicates by declared name: if your project pulls in other
subprojects that *also* `FetchContent_Declare(Splines ...)` (directly or
transitively), it's fetched, configured, and built exactly once no matter
how many times `FetchContent_MakeAvailable(Splines)` is called. Version
mismatches are resolved silently — CMake uses whichever `FetchContent_Declare`
call it processes *first* in configure order. If subprojects might reasonably
disagree on the version, declare `Splines` in your own top-level
`CMakeLists.txt` **before** `add_subdirectory`/`FetchContent`-ing anything
that also depends on it, so your pin is the one that deterministically wins.

### Installing without CMake integration

`cmake --install` (see [Building](#building)) places the compiled library
under `lib/` and the public headers under `include/Splines/` in the given
prefix, for consumption from a non-CMake build (or a manual
`-I`/`-L`/`-l` setup). There's deliberately no `SplinesConfig.cmake` /
`find_package(Splines)` support: `Splines::Splines` publicly links
`GenericContainer`, `UtilsLite`, and `quarticRootsFlocke` targets that are
built in-tree via `FetchContent` and don't export themselves when nested as
someone else's subdependency, so a CMake package config for Splines can't
resolve them either. `FetchContent` (above) is the supported way to consume
Splines from another CMake project.

## The C interface

`Splines/SplinesCinterface.h` exposes a handle-based C API
(`Spline_new`, `Spline_build`, `Spline_eval`, ...) for embedding in
non-C++ hosts — it's part of the same `Splines::Splines` target, no
separate library to link. The C API is internally synchronized (each call
is atomic), so it is memory-safe under concurrent use; it remains logically
stateful, so heavy concurrent evaluation should go through the C++ API.

## Performance

Evaluating an array of abscissas through the batch API
(`eval`/`D`/`DD( span, span )`) is **3–4× faster** than a per-point scalar
loop when the abscissas are **sorted** (resampling, plotting, arc-length
reparametrization): the batch reuses the current segment instead of
re-searching for every point. On unsorted input it auto-detects the ordering
and falls back to a plain per-point loop, staying at scalar speed — never a
regression. Batch results are always bit-for-bit identical to the scalar path.

Sorted queries, Apple M1 Max, 1000 knots, 2M queries, ns per evaluation:

| Spline  | scalar loop | batch | speedup |
|---------|------------:|------:|--------:|
| Linear  |  8.9 | 2.3 | 3.8× |
| Cubic   | 10.5 | 2.6 | 4.0× |
| Akima   | 10.4 | 2.6 | 4.0× |
| PCHIP   | 10.4 | 2.6 | 4.0× |
| Quintic | 12.5 | 4.1 | 3.0× |

Numbers are machine-dependent; reproduce with
`./bin/bench_eval [knots] [queries] [reps]` (built when
`SPLINES_BUILD_BENCHMARKS` is on).

## Testing

```sh
ctest --test-dir build --output-on-failure
```

Or just:

```sh
./build.sh test
```

`test00` through `test18` cover each curve/surface type, the C interface,
and the `nlohmann::json` / `GenericContainer` bridges (`test18`).

## License

See [license.txt](license.txt).
