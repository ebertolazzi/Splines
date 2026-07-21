# Splines MATLAB toolbox

From the repository root, configure the main project first. This resolves each
dependency from a sibling checkout or FetchContent and populates `toolbox/src`:

```sh
cmake -S . -B build/toolbox-dependencies \
  -DSPLINES_POPULATE_TOOLBOX=ON \
  -DSPLINES_UPDATE_3RDPARTY=OFF \
  -DSPLINES_INSTALL=OFF \
  -DBUILD_TESTING=OFF
```

Then configure and build the MEX targets:

```sh
cmake -S toolbox -B toolbox/build
cmake --build toolbox/build --parallel
```

As a convenience, `ruby toolbox/build.rb` runs both CMake stages and starts
from a clean `toolbox/build` directory.

To check the library in MATLAB run

```
splines_setup
```

before using it.
