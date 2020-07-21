# Python Wrapper

The python wrapper is one to one with the Splines library. 

## Requirements

The wrapper uses the Github version of `pybind11`, and it is not compatible with 
the `pip` version (that lacks cmake support).

In order to compile `pybind11` (commands valid on a Unix machine, Win version should
apply some minor modifications):

``` bash
export PYTHON_TARGET=/usr/bin/python3 # or /usr/bin/python2 for py2 version
git clone https://github.com/pybind/pybind11.git
cd pybind11
mkdir build && cd build
cmake -DPYTHON_EXECUTABLE:FILEPATH=$PYTHON_TARGET ..
sudo make install -j8
```

## Compile the wrapper

It is possible to compile the wrapper for both python 2 and 3. In order to compile
the wrapper:

``` bash
export PYTHON_TARGET=/usr/bin/python3 # or /usr/bin/python2 for py2 version
cd src_py
mkdir build && cd build
cmake -DPYTHON_EXECUTABLE:FILEPATH=$PYTHON_TARGET ..
make install -j8
```

The shared object to be imported in python is stored in `{splines_dir}/lib/lib/Splines.{ext}`,
near the Splines static library.

## Usage

In python:

``` python
import sys
sys.path.insert(0, "{splines_dir}/lib/lib")

import Splines

help(Splines)
```

The wrapper wraps one to one the original C++ library, the same interface should be expected.