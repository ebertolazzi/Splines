# get the type of OS currently running
OS=$(shell uname)
PWD=$(shell pwd)

LIB_SPLINE = libSplines.a

CC   = gcc
CXX  = g++ -std=c++11
INC  = -Isrc -Iinclude \
       -Isubmodules/GenericContainer/src \
			 -Isubmodules/Utils/src \
			 -Isubmodules/Utils/src/Utils \
			 -Isubmodules/quarticRootsFlocke/src
LIBS = -Llib -Lsubmodules/GenericContainer/lib -lSplines
DEFS =
MAKE = make

# check if the OS string contains 'Linux'
ifneq (,$(findstring Linux, $(OS)))
  WARN = -Wall
  # activate C++11 for g++ >= 4.9
  VERSION  = $(shell $(CC) -dumpversion)
  CC      += $(WARN)
  CXX     += $(WARN)
  LIBS    += -ldl
  CXXFLAGS = -pthread -Wall -O2 -fPIC -Wno-sign-compare
  AR       = ar rcs
endif

# check if the OS string contains 'Darwin'
ifneq (,$(findstring Darwin, $(OS)))
  WARN     = -Wall
  CC       = clang
  CXX      = clang++ -std=c++11 -stdlib=libc++
  VERSION  = $(shell $(CC) --version 2>&1 | grep -o "Apple LLVM version [0-9]\.[0-9]\.[0-9]" | grep -o " [0-9]\.")
  CC      += $(WARN)
  CXX     += $(WARN)
  CXXFLAGS = -Wall -O2 -fPIC -Wno-sign-compare
  AR       = libtool -static -o
endif

SRCS  = $(shell echo src/*.cc) \
        $(shell echo submodules/Utils/src/*.cc) \
        $(shell echo submodules/Utils/src/Utils/fmt/*.cc) \
        $(shell echo submodules/quarticRootsFlocke/src/*.cc) \
        $(shell echo submodules/GenericContainer/src/*.cc)
OBJS  = $(SRCS:.cc=.o)
DEPS  = $(shell echo src/*.hh) $(shell echo submodules/Utils/src/*.h*)
MKDIR = mkdir -p

# prefix for installation, use make PREFIX=/new/prefix install
# to override
PREFIX    = /usr/local
FRAMEWORK = Splines

all: lib tests

lib: lib/$(LIB_SPLINE)

tests: lib
	@echo "compile binaries\n\n"
	mkdir -p bin
	$(CXX) $(INC) $(CXXFLAGS) -o bin/test1  tests/test1.cc  $(LIBS)
	$(CXX) $(INC) $(CXXFLAGS) -o bin/test2  tests/test2.cc  $(LIBS)
	$(CXX) $(INC) $(CXXFLAGS) -o bin/test3  tests/test3.cc  $(LIBS)
	$(CXX) $(INC) $(CXXFLAGS) -o bin/test4  tests/test4.cc  $(LIBS)
	$(CXX) $(INC) $(CXXFLAGS) -o bin/test5  tests/test5.cc  $(LIBS)
	$(CXX) $(INC) $(CXXFLAGS) -o bin/test6  tests/test6.cc  $(LIBS)
	$(CXX) $(INC) $(CXXFLAGS) -o bin/test8  tests/test8.cc  $(LIBS)
	$(CXX) $(INC) $(CXXFLAGS) -o bin/test9  tests/test9.cc  $(LIBS)
	$(CXX) $(INC) $(CXXFLAGS) -o bin/test10 tests/test10.cc $(LIBS)

travis: lib tests run

include_local:
	@rm -rf lib/include
	$(MKDIR) lib
	$(MKDIR) lib/include
	@cp -f src/*.h* lib/include

.cc.o:
	$(CXX) $(INC) $(CXXFLAGS) $(DEFS) -c $< -o $@

lib/libSplines.a: $(OBJS) include_local
	$(AR) lib/libSplines.a $(OBJS)

lib/libSplines.dylib: $(OBJS) include_local
	$(CXX) -shared -o lib/libSplines.dylib $(OBJS)

lib/libSplines.so: $(OBJS) include_local
	$(CXX) -shared -o lib/libSplines.so $(OBJS)

install_local: lib
	$(MKDIR) ./lib/include
	cp GC/lib/include/*        ./lib/include
	cp src/Splines.hh          ./lib/include
	cp src/*.hxx               ./lib/include
	cp src/SplinesCinterface.h ./lib/include

install: lib
	cp src/Splines.hh          $(PREFIX)/include
	cp src/SplinesCinterface.h $(PREFIX)/include
	cp src/*.hxx               $(PREFIX)/include
	cp lib/$(LIB_SPLINE)       $(PREFIX)/lib

install_as_framework: lib
	$(MKDIR) $(PREFIX)/include/$(FRAMEWORK)
	cp src/Splines.hh          $(PREFIX)/include/$(FRAMEWORK)
	cp src/SplinesCinterface.h $(PREFIX)/include/$(FRAMEWORK)
	cp src/*.hxx               $(PREFIX)/include/$(FRAMEWORK)
	cp lib/$(LIB_SPLINE)       $(PREFIX)/lib

run:
	@echo "run binaries\n\n"
	./bin/test1
	./bin/test2
	./bin/test3
	./bin/test4
	./bin/test5
	./bin/test6
	./bin/test8
	./bin/test9
	./bin/test10

doc:
	@cd docs_build; make; cd ..

clean:
	rm -rf lib/libSplines.* lib/libGenericContainer.* lib/include src/*.o
	rm -rf bin

depend:
	makedepend -- $(INC) $(CXXFLAGS) $(DEFS) -- $(SRCS)
# DO NOT DELETE

src/SplineAkima.o: src/Splines.hh src/SplinesConfig.hh
src/SplineAkima.o: submodules/Utils/src/Utils.hh
src/SplineAkima.o: submodules/Utils/src/Utils/Utils.hxx
src/SplineAkima.o: submodules/Utils/src/Utils/rang.hpp
src/SplineAkima.o: submodules/Utils/src/Utils/fmt/printf.h
src/SplineAkima.o: submodules/Utils/src/Utils/fmt/ostream.h
src/SplineAkima.o: submodules/Utils/src/Utils/fmt/format.h
src/SplineAkima.o: submodules/Utils/src/Utils/fmt/core.h
src/SplineAkima.o: submodules/Utils/src/Utils/fmt/chrono.h
src/SplineAkima.o: submodules/Utils/src/Utils/fmt/locale.h
src/SplineAkima.o: submodules/Utils/src/Utils/fmt/ostream.h
src/SplineAkima.o: submodules/Utils/src/Utils/zstream/izstream.hpp
src/SplineAkima.o: submodules/Utils/src/Utils/zstream/zstream_common.hpp
src/SplineAkima.o: submodules/Utils/src/Utils/zstream/izstream_impl.hpp
src/SplineAkima.o: submodules/Utils/src/Utils/zstream/izstream.hpp
src/SplineAkima.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
src/SplineAkima.o: submodules/Utils/src/Utils/zstream/ozstream_impl.hpp
src/SplineAkima.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
src/SplineAkima.o: submodules/Utils/src/Utils/Trace.hxx
src/SplineAkima.o: submodules/Utils/src/Utils/Console.hxx
src/SplineAkima.o: submodules/Utils/src/Utils/Malloc.hxx
src/SplineAkima.o: submodules/Utils/src/Utils/Numbers.hxx
src/SplineAkima.o: submodules/Utils/src/Utils/TicToc.hxx
src/SplineAkima.o: submodules/Utils/src/Utils/ThreadPool.hxx
src/SplineAkima.o: submodules/GenericContainer/src/GenericContainer.hh
src/SplineAkima.o: submodules/GenericContainer/src/GenericContainerConfig.hh
src/SplineAkima2D.o: src/Splines.hh src/SplinesConfig.hh
src/SplineAkima2D.o: submodules/Utils/src/Utils.hh
src/SplineAkima2D.o: submodules/Utils/src/Utils/Utils.hxx
src/SplineAkima2D.o: submodules/Utils/src/Utils/rang.hpp
src/SplineAkima2D.o: submodules/Utils/src/Utils/fmt/printf.h
src/SplineAkima2D.o: submodules/Utils/src/Utils/fmt/ostream.h
src/SplineAkima2D.o: submodules/Utils/src/Utils/fmt/format.h
src/SplineAkima2D.o: submodules/Utils/src/Utils/fmt/core.h
src/SplineAkima2D.o: submodules/Utils/src/Utils/fmt/chrono.h
src/SplineAkima2D.o: submodules/Utils/src/Utils/fmt/locale.h
src/SplineAkima2D.o: submodules/Utils/src/Utils/fmt/ostream.h
src/SplineAkima2D.o: submodules/Utils/src/Utils/zstream/izstream.hpp
src/SplineAkima2D.o: submodules/Utils/src/Utils/zstream/zstream_common.hpp
src/SplineAkima2D.o: submodules/Utils/src/Utils/zstream/izstream_impl.hpp
src/SplineAkima2D.o: submodules/Utils/src/Utils/zstream/izstream.hpp
src/SplineAkima2D.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
src/SplineAkima2D.o: submodules/Utils/src/Utils/zstream/ozstream_impl.hpp
src/SplineAkima2D.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
src/SplineAkima2D.o: submodules/Utils/src/Utils/Trace.hxx
src/SplineAkima2D.o: submodules/Utils/src/Utils/Console.hxx
src/SplineAkima2D.o: submodules/Utils/src/Utils/Malloc.hxx
src/SplineAkima2D.o: submodules/Utils/src/Utils/Numbers.hxx
src/SplineAkima2D.o: submodules/Utils/src/Utils/TicToc.hxx
src/SplineAkima2D.o: submodules/Utils/src/Utils/ThreadPool.hxx
src/SplineAkima2D.o: submodules/GenericContainer/src/GenericContainer.hh
src/SplineAkima2D.o: submodules/GenericContainer/src/GenericContainerConfig.hh
src/SplineBessel.o: src/Splines.hh src/SplinesConfig.hh
src/SplineBessel.o: submodules/Utils/src/Utils.hh
src/SplineBessel.o: submodules/Utils/src/Utils/Utils.hxx
src/SplineBessel.o: submodules/Utils/src/Utils/rang.hpp
src/SplineBessel.o: submodules/Utils/src/Utils/fmt/printf.h
src/SplineBessel.o: submodules/Utils/src/Utils/fmt/ostream.h
src/SplineBessel.o: submodules/Utils/src/Utils/fmt/format.h
src/SplineBessel.o: submodules/Utils/src/Utils/fmt/core.h
src/SplineBessel.o: submodules/Utils/src/Utils/fmt/chrono.h
src/SplineBessel.o: submodules/Utils/src/Utils/fmt/locale.h
src/SplineBessel.o: submodules/Utils/src/Utils/fmt/ostream.h
src/SplineBessel.o: submodules/Utils/src/Utils/zstream/izstream.hpp
src/SplineBessel.o: submodules/Utils/src/Utils/zstream/zstream_common.hpp
src/SplineBessel.o: submodules/Utils/src/Utils/zstream/izstream_impl.hpp
src/SplineBessel.o: submodules/Utils/src/Utils/zstream/izstream.hpp
src/SplineBessel.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
src/SplineBessel.o: submodules/Utils/src/Utils/zstream/ozstream_impl.hpp
src/SplineBessel.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
src/SplineBessel.o: submodules/Utils/src/Utils/Trace.hxx
src/SplineBessel.o: submodules/Utils/src/Utils/Console.hxx
src/SplineBessel.o: submodules/Utils/src/Utils/Malloc.hxx
src/SplineBessel.o: submodules/Utils/src/Utils/Numbers.hxx
src/SplineBessel.o: submodules/Utils/src/Utils/TicToc.hxx
src/SplineBessel.o: submodules/Utils/src/Utils/ThreadPool.hxx
src/SplineBessel.o: submodules/GenericContainer/src/GenericContainer.hh
src/SplineBessel.o: submodules/GenericContainer/src/GenericContainerConfig.hh
src/SplineBiCubic.o: src/Splines.hh src/SplinesConfig.hh
src/SplineBiCubic.o: submodules/Utils/src/Utils.hh
src/SplineBiCubic.o: submodules/Utils/src/Utils/Utils.hxx
src/SplineBiCubic.o: submodules/Utils/src/Utils/rang.hpp
src/SplineBiCubic.o: submodules/Utils/src/Utils/fmt/printf.h
src/SplineBiCubic.o: submodules/Utils/src/Utils/fmt/ostream.h
src/SplineBiCubic.o: submodules/Utils/src/Utils/fmt/format.h
src/SplineBiCubic.o: submodules/Utils/src/Utils/fmt/core.h
src/SplineBiCubic.o: submodules/Utils/src/Utils/fmt/chrono.h
src/SplineBiCubic.o: submodules/Utils/src/Utils/fmt/locale.h
src/SplineBiCubic.o: submodules/Utils/src/Utils/fmt/ostream.h
src/SplineBiCubic.o: submodules/Utils/src/Utils/zstream/izstream.hpp
src/SplineBiCubic.o: submodules/Utils/src/Utils/zstream/zstream_common.hpp
src/SplineBiCubic.o: submodules/Utils/src/Utils/zstream/izstream_impl.hpp
src/SplineBiCubic.o: submodules/Utils/src/Utils/zstream/izstream.hpp
src/SplineBiCubic.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
src/SplineBiCubic.o: submodules/Utils/src/Utils/zstream/ozstream_impl.hpp
src/SplineBiCubic.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
src/SplineBiCubic.o: submodules/Utils/src/Utils/Trace.hxx
src/SplineBiCubic.o: submodules/Utils/src/Utils/Console.hxx
src/SplineBiCubic.o: submodules/Utils/src/Utils/Malloc.hxx
src/SplineBiCubic.o: submodules/Utils/src/Utils/Numbers.hxx
src/SplineBiCubic.o: submodules/Utils/src/Utils/TicToc.hxx
src/SplineBiCubic.o: submodules/Utils/src/Utils/ThreadPool.hxx
src/SplineBiCubic.o: submodules/GenericContainer/src/GenericContainer.hh
src/SplineBiCubic.o: submodules/GenericContainer/src/GenericContainerConfig.hh
src/SplineBiQuintic.o: src/Splines.hh src/SplinesConfig.hh
src/SplineBiQuintic.o: submodules/Utils/src/Utils.hh
src/SplineBiQuintic.o: submodules/Utils/src/Utils/Utils.hxx
src/SplineBiQuintic.o: submodules/Utils/src/Utils/rang.hpp
src/SplineBiQuintic.o: submodules/Utils/src/Utils/fmt/printf.h
src/SplineBiQuintic.o: submodules/Utils/src/Utils/fmt/ostream.h
src/SplineBiQuintic.o: submodules/Utils/src/Utils/fmt/format.h
src/SplineBiQuintic.o: submodules/Utils/src/Utils/fmt/core.h
src/SplineBiQuintic.o: submodules/Utils/src/Utils/fmt/chrono.h
src/SplineBiQuintic.o: submodules/Utils/src/Utils/fmt/locale.h
src/SplineBiQuintic.o: submodules/Utils/src/Utils/fmt/ostream.h
src/SplineBiQuintic.o: submodules/Utils/src/Utils/zstream/izstream.hpp
src/SplineBiQuintic.o: submodules/Utils/src/Utils/zstream/zstream_common.hpp
src/SplineBiQuintic.o: submodules/Utils/src/Utils/zstream/izstream_impl.hpp
src/SplineBiQuintic.o: submodules/Utils/src/Utils/zstream/izstream.hpp
src/SplineBiQuintic.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
src/SplineBiQuintic.o: submodules/Utils/src/Utils/zstream/ozstream_impl.hpp
src/SplineBiQuintic.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
src/SplineBiQuintic.o: submodules/Utils/src/Utils/Trace.hxx
src/SplineBiQuintic.o: submodules/Utils/src/Utils/Console.hxx
src/SplineBiQuintic.o: submodules/Utils/src/Utils/Malloc.hxx
src/SplineBiQuintic.o: submodules/Utils/src/Utils/Numbers.hxx
src/SplineBiQuintic.o: submodules/Utils/src/Utils/TicToc.hxx
src/SplineBiQuintic.o: submodules/Utils/src/Utils/ThreadPool.hxx
src/SplineBiQuintic.o: submodules/GenericContainer/src/GenericContainer.hh
src/SplineBiQuintic.o: submodules/GenericContainer/src/GenericContainerConfig.hh
src/SplineBilinear.o: src/Splines.hh src/SplinesConfig.hh
src/SplineBilinear.o: submodules/Utils/src/Utils.hh
src/SplineBilinear.o: submodules/Utils/src/Utils/Utils.hxx
src/SplineBilinear.o: submodules/Utils/src/Utils/rang.hpp
src/SplineBilinear.o: submodules/Utils/src/Utils/fmt/printf.h
src/SplineBilinear.o: submodules/Utils/src/Utils/fmt/ostream.h
src/SplineBilinear.o: submodules/Utils/src/Utils/fmt/format.h
src/SplineBilinear.o: submodules/Utils/src/Utils/fmt/core.h
src/SplineBilinear.o: submodules/Utils/src/Utils/fmt/chrono.h
src/SplineBilinear.o: submodules/Utils/src/Utils/fmt/locale.h
src/SplineBilinear.o: submodules/Utils/src/Utils/fmt/ostream.h
src/SplineBilinear.o: submodules/Utils/src/Utils/zstream/izstream.hpp
src/SplineBilinear.o: submodules/Utils/src/Utils/zstream/zstream_common.hpp
src/SplineBilinear.o: submodules/Utils/src/Utils/zstream/izstream_impl.hpp
src/SplineBilinear.o: submodules/Utils/src/Utils/zstream/izstream.hpp
src/SplineBilinear.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
src/SplineBilinear.o: submodules/Utils/src/Utils/zstream/ozstream_impl.hpp
src/SplineBilinear.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
src/SplineBilinear.o: submodules/Utils/src/Utils/Trace.hxx
src/SplineBilinear.o: submodules/Utils/src/Utils/Console.hxx
src/SplineBilinear.o: submodules/Utils/src/Utils/Malloc.hxx
src/SplineBilinear.o: submodules/Utils/src/Utils/Numbers.hxx
src/SplineBilinear.o: submodules/Utils/src/Utils/TicToc.hxx
src/SplineBilinear.o: submodules/Utils/src/Utils/ThreadPool.hxx
src/SplineBilinear.o: submodules/GenericContainer/src/GenericContainer.hh
src/SplineBilinear.o: submodules/GenericContainer/src/GenericContainerConfig.hh
src/SplineConstant.o: src/Splines.hh src/SplinesConfig.hh
src/SplineConstant.o: submodules/Utils/src/Utils.hh
src/SplineConstant.o: submodules/Utils/src/Utils/Utils.hxx
src/SplineConstant.o: submodules/Utils/src/Utils/rang.hpp
src/SplineConstant.o: submodules/Utils/src/Utils/fmt/printf.h
src/SplineConstant.o: submodules/Utils/src/Utils/fmt/ostream.h
src/SplineConstant.o: submodules/Utils/src/Utils/fmt/format.h
src/SplineConstant.o: submodules/Utils/src/Utils/fmt/core.h
src/SplineConstant.o: submodules/Utils/src/Utils/fmt/chrono.h
src/SplineConstant.o: submodules/Utils/src/Utils/fmt/locale.h
src/SplineConstant.o: submodules/Utils/src/Utils/fmt/ostream.h
src/SplineConstant.o: submodules/Utils/src/Utils/zstream/izstream.hpp
src/SplineConstant.o: submodules/Utils/src/Utils/zstream/zstream_common.hpp
src/SplineConstant.o: submodules/Utils/src/Utils/zstream/izstream_impl.hpp
src/SplineConstant.o: submodules/Utils/src/Utils/zstream/izstream.hpp
src/SplineConstant.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
src/SplineConstant.o: submodules/Utils/src/Utils/zstream/ozstream_impl.hpp
src/SplineConstant.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
src/SplineConstant.o: submodules/Utils/src/Utils/Trace.hxx
src/SplineConstant.o: submodules/Utils/src/Utils/Console.hxx
src/SplineConstant.o: submodules/Utils/src/Utils/Malloc.hxx
src/SplineConstant.o: submodules/Utils/src/Utils/Numbers.hxx
src/SplineConstant.o: submodules/Utils/src/Utils/TicToc.hxx
src/SplineConstant.o: submodules/Utils/src/Utils/ThreadPool.hxx
src/SplineConstant.o: submodules/GenericContainer/src/GenericContainer.hh
src/SplineConstant.o: submodules/GenericContainer/src/GenericContainerConfig.hh
src/SplineCubic.o: src/Splines.hh src/SplinesConfig.hh
src/SplineCubic.o: submodules/Utils/src/Utils.hh
src/SplineCubic.o: submodules/Utils/src/Utils/Utils.hxx
src/SplineCubic.o: submodules/Utils/src/Utils/rang.hpp
src/SplineCubic.o: submodules/Utils/src/Utils/fmt/printf.h
src/SplineCubic.o: submodules/Utils/src/Utils/fmt/ostream.h
src/SplineCubic.o: submodules/Utils/src/Utils/fmt/format.h
src/SplineCubic.o: submodules/Utils/src/Utils/fmt/core.h
src/SplineCubic.o: submodules/Utils/src/Utils/fmt/chrono.h
src/SplineCubic.o: submodules/Utils/src/Utils/fmt/locale.h
src/SplineCubic.o: submodules/Utils/src/Utils/fmt/ostream.h
src/SplineCubic.o: submodules/Utils/src/Utils/zstream/izstream.hpp
src/SplineCubic.o: submodules/Utils/src/Utils/zstream/zstream_common.hpp
src/SplineCubic.o: submodules/Utils/src/Utils/zstream/izstream_impl.hpp
src/SplineCubic.o: submodules/Utils/src/Utils/zstream/izstream.hpp
src/SplineCubic.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
src/SplineCubic.o: submodules/Utils/src/Utils/zstream/ozstream_impl.hpp
src/SplineCubic.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
src/SplineCubic.o: submodules/Utils/src/Utils/Trace.hxx
src/SplineCubic.o: submodules/Utils/src/Utils/Console.hxx
src/SplineCubic.o: submodules/Utils/src/Utils/Malloc.hxx
src/SplineCubic.o: submodules/Utils/src/Utils/Numbers.hxx
src/SplineCubic.o: submodules/Utils/src/Utils/TicToc.hxx
src/SplineCubic.o: submodules/Utils/src/Utils/ThreadPool.hxx
src/SplineCubic.o: submodules/GenericContainer/src/GenericContainer.hh
src/SplineCubic.o: submodules/GenericContainer/src/GenericContainerConfig.hh
src/SplineCubicBase.o: src/Splines.hh src/SplinesConfig.hh
src/SplineCubicBase.o: submodules/Utils/src/Utils.hh
src/SplineCubicBase.o: submodules/Utils/src/Utils/Utils.hxx
src/SplineCubicBase.o: submodules/Utils/src/Utils/rang.hpp
src/SplineCubicBase.o: submodules/Utils/src/Utils/fmt/printf.h
src/SplineCubicBase.o: submodules/Utils/src/Utils/fmt/ostream.h
src/SplineCubicBase.o: submodules/Utils/src/Utils/fmt/format.h
src/SplineCubicBase.o: submodules/Utils/src/Utils/fmt/core.h
src/SplineCubicBase.o: submodules/Utils/src/Utils/fmt/chrono.h
src/SplineCubicBase.o: submodules/Utils/src/Utils/fmt/locale.h
src/SplineCubicBase.o: submodules/Utils/src/Utils/fmt/ostream.h
src/SplineCubicBase.o: submodules/Utils/src/Utils/zstream/izstream.hpp
src/SplineCubicBase.o: submodules/Utils/src/Utils/zstream/zstream_common.hpp
src/SplineCubicBase.o: submodules/Utils/src/Utils/zstream/izstream_impl.hpp
src/SplineCubicBase.o: submodules/Utils/src/Utils/zstream/izstream.hpp
src/SplineCubicBase.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
src/SplineCubicBase.o: submodules/Utils/src/Utils/zstream/ozstream_impl.hpp
src/SplineCubicBase.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
src/SplineCubicBase.o: submodules/Utils/src/Utils/Trace.hxx
src/SplineCubicBase.o: submodules/Utils/src/Utils/Console.hxx
src/SplineCubicBase.o: submodules/Utils/src/Utils/Malloc.hxx
src/SplineCubicBase.o: submodules/Utils/src/Utils/Numbers.hxx
src/SplineCubicBase.o: submodules/Utils/src/Utils/TicToc.hxx
src/SplineCubicBase.o: submodules/Utils/src/Utils/ThreadPool.hxx
src/SplineCubicBase.o: submodules/GenericContainer/src/GenericContainer.hh
src/SplineCubicBase.o: submodules/GenericContainer/src/GenericContainerConfig.hh
src/SplineHermite.o: src/Splines.hh src/SplinesConfig.hh
src/SplineHermite.o: submodules/Utils/src/Utils.hh
src/SplineHermite.o: submodules/Utils/src/Utils/Utils.hxx
src/SplineHermite.o: submodules/Utils/src/Utils/rang.hpp
src/SplineHermite.o: submodules/Utils/src/Utils/fmt/printf.h
src/SplineHermite.o: submodules/Utils/src/Utils/fmt/ostream.h
src/SplineHermite.o: submodules/Utils/src/Utils/fmt/format.h
src/SplineHermite.o: submodules/Utils/src/Utils/fmt/core.h
src/SplineHermite.o: submodules/Utils/src/Utils/fmt/chrono.h
src/SplineHermite.o: submodules/Utils/src/Utils/fmt/locale.h
src/SplineHermite.o: submodules/Utils/src/Utils/fmt/ostream.h
src/SplineHermite.o: submodules/Utils/src/Utils/zstream/izstream.hpp
src/SplineHermite.o: submodules/Utils/src/Utils/zstream/zstream_common.hpp
src/SplineHermite.o: submodules/Utils/src/Utils/zstream/izstream_impl.hpp
src/SplineHermite.o: submodules/Utils/src/Utils/zstream/izstream.hpp
src/SplineHermite.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
src/SplineHermite.o: submodules/Utils/src/Utils/zstream/ozstream_impl.hpp
src/SplineHermite.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
src/SplineHermite.o: submodules/Utils/src/Utils/Trace.hxx
src/SplineHermite.o: submodules/Utils/src/Utils/Console.hxx
src/SplineHermite.o: submodules/Utils/src/Utils/Malloc.hxx
src/SplineHermite.o: submodules/Utils/src/Utils/Numbers.hxx
src/SplineHermite.o: submodules/Utils/src/Utils/TicToc.hxx
src/SplineHermite.o: submodules/Utils/src/Utils/ThreadPool.hxx
src/SplineHermite.o: submodules/GenericContainer/src/GenericContainer.hh
src/SplineHermite.o: submodules/GenericContainer/src/GenericContainerConfig.hh
src/SplineLinear.o: src/Splines.hh src/SplinesConfig.hh
src/SplineLinear.o: submodules/Utils/src/Utils.hh
src/SplineLinear.o: submodules/Utils/src/Utils/Utils.hxx
src/SplineLinear.o: submodules/Utils/src/Utils/rang.hpp
src/SplineLinear.o: submodules/Utils/src/Utils/fmt/printf.h
src/SplineLinear.o: submodules/Utils/src/Utils/fmt/ostream.h
src/SplineLinear.o: submodules/Utils/src/Utils/fmt/format.h
src/SplineLinear.o: submodules/Utils/src/Utils/fmt/core.h
src/SplineLinear.o: submodules/Utils/src/Utils/fmt/chrono.h
src/SplineLinear.o: submodules/Utils/src/Utils/fmt/locale.h
src/SplineLinear.o: submodules/Utils/src/Utils/fmt/ostream.h
src/SplineLinear.o: submodules/Utils/src/Utils/zstream/izstream.hpp
src/SplineLinear.o: submodules/Utils/src/Utils/zstream/zstream_common.hpp
src/SplineLinear.o: submodules/Utils/src/Utils/zstream/izstream_impl.hpp
src/SplineLinear.o: submodules/Utils/src/Utils/zstream/izstream.hpp
src/SplineLinear.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
src/SplineLinear.o: submodules/Utils/src/Utils/zstream/ozstream_impl.hpp
src/SplineLinear.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
src/SplineLinear.o: submodules/Utils/src/Utils/Trace.hxx
src/SplineLinear.o: submodules/Utils/src/Utils/Console.hxx
src/SplineLinear.o: submodules/Utils/src/Utils/Malloc.hxx
src/SplineLinear.o: submodules/Utils/src/Utils/Numbers.hxx
src/SplineLinear.o: submodules/Utils/src/Utils/TicToc.hxx
src/SplineLinear.o: submodules/Utils/src/Utils/ThreadPool.hxx
src/SplineLinear.o: submodules/GenericContainer/src/GenericContainer.hh
src/SplineLinear.o: submodules/GenericContainer/src/GenericContainerConfig.hh
src/SplinePchip.o: src/Splines.hh src/SplinesConfig.hh
src/SplinePchip.o: submodules/Utils/src/Utils.hh
src/SplinePchip.o: submodules/Utils/src/Utils/Utils.hxx
src/SplinePchip.o: submodules/Utils/src/Utils/rang.hpp
src/SplinePchip.o: submodules/Utils/src/Utils/fmt/printf.h
src/SplinePchip.o: submodules/Utils/src/Utils/fmt/ostream.h
src/SplinePchip.o: submodules/Utils/src/Utils/fmt/format.h
src/SplinePchip.o: submodules/Utils/src/Utils/fmt/core.h
src/SplinePchip.o: submodules/Utils/src/Utils/fmt/chrono.h
src/SplinePchip.o: submodules/Utils/src/Utils/fmt/locale.h
src/SplinePchip.o: submodules/Utils/src/Utils/fmt/ostream.h
src/SplinePchip.o: submodules/Utils/src/Utils/zstream/izstream.hpp
src/SplinePchip.o: submodules/Utils/src/Utils/zstream/zstream_common.hpp
src/SplinePchip.o: submodules/Utils/src/Utils/zstream/izstream_impl.hpp
src/SplinePchip.o: submodules/Utils/src/Utils/zstream/izstream.hpp
src/SplinePchip.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
src/SplinePchip.o: submodules/Utils/src/Utils/zstream/ozstream_impl.hpp
src/SplinePchip.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
src/SplinePchip.o: submodules/Utils/src/Utils/Trace.hxx
src/SplinePchip.o: submodules/Utils/src/Utils/Console.hxx
src/SplinePchip.o: submodules/Utils/src/Utils/Malloc.hxx
src/SplinePchip.o: submodules/Utils/src/Utils/Numbers.hxx
src/SplinePchip.o: submodules/Utils/src/Utils/TicToc.hxx
src/SplinePchip.o: submodules/Utils/src/Utils/ThreadPool.hxx
src/SplinePchip.o: submodules/GenericContainer/src/GenericContainer.hh
src/SplinePchip.o: submodules/GenericContainer/src/GenericContainerConfig.hh
src/SplinePchip.o: src/SplinesUtils.hh
src/SplineQuintic.o: src/Splines.hh src/SplinesConfig.hh
src/SplineQuintic.o: submodules/Utils/src/Utils.hh
src/SplineQuintic.o: submodules/Utils/src/Utils/Utils.hxx
src/SplineQuintic.o: submodules/Utils/src/Utils/rang.hpp
src/SplineQuintic.o: submodules/Utils/src/Utils/fmt/printf.h
src/SplineQuintic.o: submodules/Utils/src/Utils/fmt/ostream.h
src/SplineQuintic.o: submodules/Utils/src/Utils/fmt/format.h
src/SplineQuintic.o: submodules/Utils/src/Utils/fmt/core.h
src/SplineQuintic.o: submodules/Utils/src/Utils/fmt/chrono.h
src/SplineQuintic.o: submodules/Utils/src/Utils/fmt/locale.h
src/SplineQuintic.o: submodules/Utils/src/Utils/fmt/ostream.h
src/SplineQuintic.o: submodules/Utils/src/Utils/zstream/izstream.hpp
src/SplineQuintic.o: submodules/Utils/src/Utils/zstream/zstream_common.hpp
src/SplineQuintic.o: submodules/Utils/src/Utils/zstream/izstream_impl.hpp
src/SplineQuintic.o: submodules/Utils/src/Utils/zstream/izstream.hpp
src/SplineQuintic.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
src/SplineQuintic.o: submodules/Utils/src/Utils/zstream/ozstream_impl.hpp
src/SplineQuintic.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
src/SplineQuintic.o: submodules/Utils/src/Utils/Trace.hxx
src/SplineQuintic.o: submodules/Utils/src/Utils/Console.hxx
src/SplineQuintic.o: submodules/Utils/src/Utils/Malloc.hxx
src/SplineQuintic.o: submodules/Utils/src/Utils/Numbers.hxx
src/SplineQuintic.o: submodules/Utils/src/Utils/TicToc.hxx
src/SplineQuintic.o: submodules/Utils/src/Utils/ThreadPool.hxx
src/SplineQuintic.o: submodules/GenericContainer/src/GenericContainer.hh
src/SplineQuintic.o: submodules/GenericContainer/src/GenericContainerConfig.hh
src/SplineQuintic.o: src/SplinesUtils.hh
src/SplineQuinticBase.o: src/Splines.hh src/SplinesConfig.hh
src/SplineQuinticBase.o: submodules/Utils/src/Utils.hh
src/SplineQuinticBase.o: submodules/Utils/src/Utils/Utils.hxx
src/SplineQuinticBase.o: submodules/Utils/src/Utils/rang.hpp
src/SplineQuinticBase.o: submodules/Utils/src/Utils/fmt/printf.h
src/SplineQuinticBase.o: submodules/Utils/src/Utils/fmt/ostream.h
src/SplineQuinticBase.o: submodules/Utils/src/Utils/fmt/format.h
src/SplineQuinticBase.o: submodules/Utils/src/Utils/fmt/core.h
src/SplineQuinticBase.o: submodules/Utils/src/Utils/fmt/chrono.h
src/SplineQuinticBase.o: submodules/Utils/src/Utils/fmt/locale.h
src/SplineQuinticBase.o: submodules/Utils/src/Utils/fmt/ostream.h
src/SplineQuinticBase.o: submodules/Utils/src/Utils/zstream/izstream.hpp
src/SplineQuinticBase.o: submodules/Utils/src/Utils/zstream/zstream_common.hpp
src/SplineQuinticBase.o: submodules/Utils/src/Utils/zstream/izstream_impl.hpp
src/SplineQuinticBase.o: submodules/Utils/src/Utils/zstream/izstream.hpp
src/SplineQuinticBase.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
src/SplineQuinticBase.o: submodules/Utils/src/Utils/zstream/ozstream_impl.hpp
src/SplineQuinticBase.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
src/SplineQuinticBase.o: submodules/Utils/src/Utils/Trace.hxx
src/SplineQuinticBase.o: submodules/Utils/src/Utils/Console.hxx
src/SplineQuinticBase.o: submodules/Utils/src/Utils/Malloc.hxx
src/SplineQuinticBase.o: submodules/Utils/src/Utils/Numbers.hxx
src/SplineQuinticBase.o: submodules/Utils/src/Utils/TicToc.hxx
src/SplineQuinticBase.o: submodules/Utils/src/Utils/ThreadPool.hxx
src/SplineQuinticBase.o: submodules/GenericContainer/src/GenericContainer.hh
src/SplineQuinticBase.o: submodules/GenericContainer/src/GenericContainerConfig.hh
src/SplineSet.o: src/Splines.hh src/SplinesConfig.hh
src/SplineSet.o: submodules/Utils/src/Utils.hh
src/SplineSet.o: submodules/Utils/src/Utils/Utils.hxx
src/SplineSet.o: submodules/Utils/src/Utils/rang.hpp
src/SplineSet.o: submodules/Utils/src/Utils/fmt/printf.h
src/SplineSet.o: submodules/Utils/src/Utils/fmt/ostream.h
src/SplineSet.o: submodules/Utils/src/Utils/fmt/format.h
src/SplineSet.o: submodules/Utils/src/Utils/fmt/core.h
src/SplineSet.o: submodules/Utils/src/Utils/fmt/chrono.h
src/SplineSet.o: submodules/Utils/src/Utils/fmt/locale.h
src/SplineSet.o: submodules/Utils/src/Utils/fmt/ostream.h
src/SplineSet.o: submodules/Utils/src/Utils/zstream/izstream.hpp
src/SplineSet.o: submodules/Utils/src/Utils/zstream/zstream_common.hpp
src/SplineSet.o: submodules/Utils/src/Utils/zstream/izstream_impl.hpp
src/SplineSet.o: submodules/Utils/src/Utils/zstream/izstream.hpp
src/SplineSet.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
src/SplineSet.o: submodules/Utils/src/Utils/zstream/ozstream_impl.hpp
src/SplineSet.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
src/SplineSet.o: submodules/Utils/src/Utils/Trace.hxx
src/SplineSet.o: submodules/Utils/src/Utils/Console.hxx
src/SplineSet.o: submodules/Utils/src/Utils/Malloc.hxx
src/SplineSet.o: submodules/Utils/src/Utils/Numbers.hxx
src/SplineSet.o: submodules/Utils/src/Utils/TicToc.hxx
src/SplineSet.o: submodules/Utils/src/Utils/ThreadPool.hxx
src/SplineSet.o: submodules/GenericContainer/src/GenericContainer.hh
src/SplineSet.o: submodules/GenericContainer/src/GenericContainerConfig.hh
src/SplineSet.o: src/SplinesUtils.hh
src/SplineSet.o: submodules/quarticRootsFlocke/src/PolynomialRoots.hh
src/SplineSetGC.o: src/Splines.hh src/SplinesConfig.hh
src/SplineSetGC.o: submodules/Utils/src/Utils.hh
src/SplineSetGC.o: submodules/Utils/src/Utils/Utils.hxx
src/SplineSetGC.o: submodules/Utils/src/Utils/rang.hpp
src/SplineSetGC.o: submodules/Utils/src/Utils/fmt/printf.h
src/SplineSetGC.o: submodules/Utils/src/Utils/fmt/ostream.h
src/SplineSetGC.o: submodules/Utils/src/Utils/fmt/format.h
src/SplineSetGC.o: submodules/Utils/src/Utils/fmt/core.h
src/SplineSetGC.o: submodules/Utils/src/Utils/fmt/chrono.h
src/SplineSetGC.o: submodules/Utils/src/Utils/fmt/locale.h
src/SplineSetGC.o: submodules/Utils/src/Utils/fmt/ostream.h
src/SplineSetGC.o: submodules/Utils/src/Utils/zstream/izstream.hpp
src/SplineSetGC.o: submodules/Utils/src/Utils/zstream/zstream_common.hpp
src/SplineSetGC.o: submodules/Utils/src/Utils/zstream/izstream_impl.hpp
src/SplineSetGC.o: submodules/Utils/src/Utils/zstream/izstream.hpp
src/SplineSetGC.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
src/SplineSetGC.o: submodules/Utils/src/Utils/zstream/ozstream_impl.hpp
src/SplineSetGC.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
src/SplineSetGC.o: submodules/Utils/src/Utils/Trace.hxx
src/SplineSetGC.o: submodules/Utils/src/Utils/Console.hxx
src/SplineSetGC.o: submodules/Utils/src/Utils/Malloc.hxx
src/SplineSetGC.o: submodules/Utils/src/Utils/Numbers.hxx
src/SplineSetGC.o: submodules/Utils/src/Utils/TicToc.hxx
src/SplineSetGC.o: submodules/Utils/src/Utils/ThreadPool.hxx
src/SplineSetGC.o: submodules/GenericContainer/src/GenericContainer.hh
src/SplineSetGC.o: submodules/GenericContainer/src/GenericContainerConfig.hh
src/SplineVec.o: src/Splines.hh src/SplinesConfig.hh
src/SplineVec.o: submodules/Utils/src/Utils.hh
src/SplineVec.o: submodules/Utils/src/Utils/Utils.hxx
src/SplineVec.o: submodules/Utils/src/Utils/rang.hpp
src/SplineVec.o: submodules/Utils/src/Utils/fmt/printf.h
src/SplineVec.o: submodules/Utils/src/Utils/fmt/ostream.h
src/SplineVec.o: submodules/Utils/src/Utils/fmt/format.h
src/SplineVec.o: submodules/Utils/src/Utils/fmt/core.h
src/SplineVec.o: submodules/Utils/src/Utils/fmt/chrono.h
src/SplineVec.o: submodules/Utils/src/Utils/fmt/locale.h
src/SplineVec.o: submodules/Utils/src/Utils/fmt/ostream.h
src/SplineVec.o: submodules/Utils/src/Utils/zstream/izstream.hpp
src/SplineVec.o: submodules/Utils/src/Utils/zstream/zstream_common.hpp
src/SplineVec.o: submodules/Utils/src/Utils/zstream/izstream_impl.hpp
src/SplineVec.o: submodules/Utils/src/Utils/zstream/izstream.hpp
src/SplineVec.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
src/SplineVec.o: submodules/Utils/src/Utils/zstream/ozstream_impl.hpp
src/SplineVec.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
src/SplineVec.o: submodules/Utils/src/Utils/Trace.hxx
src/SplineVec.o: submodules/Utils/src/Utils/Console.hxx
src/SplineVec.o: submodules/Utils/src/Utils/Malloc.hxx
src/SplineVec.o: submodules/Utils/src/Utils/Numbers.hxx
src/SplineVec.o: submodules/Utils/src/Utils/TicToc.hxx
src/SplineVec.o: submodules/Utils/src/Utils/ThreadPool.hxx
src/SplineVec.o: submodules/GenericContainer/src/GenericContainer.hh
src/SplineVec.o: submodules/GenericContainer/src/GenericContainerConfig.hh
src/Splines.o: src/SplinesUtils.hh src/Splines.hh src/SplinesConfig.hh
src/Splines.o: submodules/Utils/src/Utils.hh
src/Splines.o: submodules/Utils/src/Utils/Utils.hxx
src/Splines.o: submodules/Utils/src/Utils/rang.hpp
src/Splines.o: submodules/Utils/src/Utils/fmt/printf.h
src/Splines.o: submodules/Utils/src/Utils/fmt/ostream.h
src/Splines.o: submodules/Utils/src/Utils/fmt/format.h
src/Splines.o: submodules/Utils/src/Utils/fmt/core.h
src/Splines.o: submodules/Utils/src/Utils/fmt/chrono.h
src/Splines.o: submodules/Utils/src/Utils/fmt/locale.h
src/Splines.o: submodules/Utils/src/Utils/fmt/ostream.h
src/Splines.o: submodules/Utils/src/Utils/zstream/izstream.hpp
src/Splines.o: submodules/Utils/src/Utils/zstream/zstream_common.hpp
src/Splines.o: submodules/Utils/src/Utils/zstream/izstream_impl.hpp
src/Splines.o: submodules/Utils/src/Utils/zstream/izstream.hpp
src/Splines.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
src/Splines.o: submodules/Utils/src/Utils/zstream/ozstream_impl.hpp
src/Splines.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
src/Splines.o: submodules/Utils/src/Utils/Trace.hxx
src/Splines.o: submodules/Utils/src/Utils/Console.hxx
src/Splines.o: submodules/Utils/src/Utils/Malloc.hxx
src/Splines.o: submodules/Utils/src/Utils/Numbers.hxx
src/Splines.o: submodules/Utils/src/Utils/TicToc.hxx
src/Splines.o: submodules/Utils/src/Utils/ThreadPool.hxx
src/Splines.o: submodules/GenericContainer/src/GenericContainer.hh
src/Splines.o: submodules/GenericContainer/src/GenericContainerConfig.hh
src/Splines1D.o: src/Splines.hh src/SplinesConfig.hh
src/Splines1D.o: submodules/Utils/src/Utils.hh
src/Splines1D.o: submodules/Utils/src/Utils/Utils.hxx
src/Splines1D.o: submodules/Utils/src/Utils/rang.hpp
src/Splines1D.o: submodules/Utils/src/Utils/fmt/printf.h
src/Splines1D.o: submodules/Utils/src/Utils/fmt/ostream.h
src/Splines1D.o: submodules/Utils/src/Utils/fmt/format.h
src/Splines1D.o: submodules/Utils/src/Utils/fmt/core.h
src/Splines1D.o: submodules/Utils/src/Utils/fmt/chrono.h
src/Splines1D.o: submodules/Utils/src/Utils/fmt/locale.h
src/Splines1D.o: submodules/Utils/src/Utils/fmt/ostream.h
src/Splines1D.o: submodules/Utils/src/Utils/zstream/izstream.hpp
src/Splines1D.o: submodules/Utils/src/Utils/zstream/zstream_common.hpp
src/Splines1D.o: submodules/Utils/src/Utils/zstream/izstream_impl.hpp
src/Splines1D.o: submodules/Utils/src/Utils/zstream/izstream.hpp
src/Splines1D.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
src/Splines1D.o: submodules/Utils/src/Utils/zstream/ozstream_impl.hpp
src/Splines1D.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
src/Splines1D.o: submodules/Utils/src/Utils/Trace.hxx
src/Splines1D.o: submodules/Utils/src/Utils/Console.hxx
src/Splines1D.o: submodules/Utils/src/Utils/Malloc.hxx
src/Splines1D.o: submodules/Utils/src/Utils/Numbers.hxx
src/Splines1D.o: submodules/Utils/src/Utils/TicToc.hxx
src/Splines1D.o: submodules/Utils/src/Utils/ThreadPool.hxx
src/Splines1D.o: submodules/GenericContainer/src/GenericContainer.hh
src/Splines1D.o: submodules/GenericContainer/src/GenericContainerConfig.hh
src/Splines2D.o: src/Splines.hh src/SplinesConfig.hh
src/Splines2D.o: submodules/Utils/src/Utils.hh
src/Splines2D.o: submodules/Utils/src/Utils/Utils.hxx
src/Splines2D.o: submodules/Utils/src/Utils/rang.hpp
src/Splines2D.o: submodules/Utils/src/Utils/fmt/printf.h
src/Splines2D.o: submodules/Utils/src/Utils/fmt/ostream.h
src/Splines2D.o: submodules/Utils/src/Utils/fmt/format.h
src/Splines2D.o: submodules/Utils/src/Utils/fmt/core.h
src/Splines2D.o: submodules/Utils/src/Utils/fmt/chrono.h
src/Splines2D.o: submodules/Utils/src/Utils/fmt/locale.h
src/Splines2D.o: submodules/Utils/src/Utils/fmt/ostream.h
src/Splines2D.o: submodules/Utils/src/Utils/zstream/izstream.hpp
src/Splines2D.o: submodules/Utils/src/Utils/zstream/zstream_common.hpp
src/Splines2D.o: submodules/Utils/src/Utils/zstream/izstream_impl.hpp
src/Splines2D.o: submodules/Utils/src/Utils/zstream/izstream.hpp
src/Splines2D.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
src/Splines2D.o: submodules/Utils/src/Utils/zstream/ozstream_impl.hpp
src/Splines2D.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
src/Splines2D.o: submodules/Utils/src/Utils/Trace.hxx
src/Splines2D.o: submodules/Utils/src/Utils/Console.hxx
src/Splines2D.o: submodules/Utils/src/Utils/Malloc.hxx
src/Splines2D.o: submodules/Utils/src/Utils/Numbers.hxx
src/Splines2D.o: submodules/Utils/src/Utils/TicToc.hxx
src/Splines2D.o: submodules/Utils/src/Utils/ThreadPool.hxx
src/Splines2D.o: submodules/GenericContainer/src/GenericContainer.hh
src/Splines2D.o: submodules/GenericContainer/src/GenericContainerConfig.hh
src/SplinesBivariate.o: src/Splines.hh src/SplinesConfig.hh
src/SplinesBivariate.o: submodules/Utils/src/Utils.hh
src/SplinesBivariate.o: submodules/Utils/src/Utils/Utils.hxx
src/SplinesBivariate.o: submodules/Utils/src/Utils/rang.hpp
src/SplinesBivariate.o: submodules/Utils/src/Utils/fmt/printf.h
src/SplinesBivariate.o: submodules/Utils/src/Utils/fmt/ostream.h
src/SplinesBivariate.o: submodules/Utils/src/Utils/fmt/format.h
src/SplinesBivariate.o: submodules/Utils/src/Utils/fmt/core.h
src/SplinesBivariate.o: submodules/Utils/src/Utils/fmt/chrono.h
src/SplinesBivariate.o: submodules/Utils/src/Utils/fmt/locale.h
src/SplinesBivariate.o: submodules/Utils/src/Utils/fmt/ostream.h
src/SplinesBivariate.o: submodules/Utils/src/Utils/zstream/izstream.hpp
src/SplinesBivariate.o: submodules/Utils/src/Utils/zstream/zstream_common.hpp
src/SplinesBivariate.o: submodules/Utils/src/Utils/zstream/izstream_impl.hpp
src/SplinesBivariate.o: submodules/Utils/src/Utils/zstream/izstream.hpp
src/SplinesBivariate.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
src/SplinesBivariate.o: submodules/Utils/src/Utils/zstream/ozstream_impl.hpp
src/SplinesBivariate.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
src/SplinesBivariate.o: submodules/Utils/src/Utils/Trace.hxx
src/SplinesBivariate.o: submodules/Utils/src/Utils/Console.hxx
src/SplinesBivariate.o: submodules/Utils/src/Utils/Malloc.hxx
src/SplinesBivariate.o: submodules/Utils/src/Utils/Numbers.hxx
src/SplinesBivariate.o: submodules/Utils/src/Utils/TicToc.hxx
src/SplinesBivariate.o: submodules/Utils/src/Utils/ThreadPool.hxx
src/SplinesBivariate.o: submodules/GenericContainer/src/GenericContainer.hh
src/SplinesBivariate.o: submodules/GenericContainer/src/GenericContainerConfig.hh
src/SplinesCinterface.o: src/Splines.hh src/SplinesConfig.hh
src/SplinesCinterface.o: submodules/Utils/src/Utils.hh
src/SplinesCinterface.o: submodules/Utils/src/Utils/Utils.hxx
src/SplinesCinterface.o: submodules/Utils/src/Utils/rang.hpp
src/SplinesCinterface.o: submodules/Utils/src/Utils/fmt/printf.h
src/SplinesCinterface.o: submodules/Utils/src/Utils/fmt/ostream.h
src/SplinesCinterface.o: submodules/Utils/src/Utils/fmt/format.h
src/SplinesCinterface.o: submodules/Utils/src/Utils/fmt/core.h
src/SplinesCinterface.o: submodules/Utils/src/Utils/fmt/chrono.h
src/SplinesCinterface.o: submodules/Utils/src/Utils/fmt/locale.h
src/SplinesCinterface.o: submodules/Utils/src/Utils/fmt/ostream.h
src/SplinesCinterface.o: submodules/Utils/src/Utils/zstream/izstream.hpp
src/SplinesCinterface.o: submodules/Utils/src/Utils/zstream/zstream_common.hpp
src/SplinesCinterface.o: submodules/Utils/src/Utils/zstream/izstream_impl.hpp
src/SplinesCinterface.o: submodules/Utils/src/Utils/zstream/izstream.hpp
src/SplinesCinterface.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
src/SplinesCinterface.o: submodules/Utils/src/Utils/zstream/ozstream_impl.hpp
src/SplinesCinterface.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
src/SplinesCinterface.o: submodules/Utils/src/Utils/Trace.hxx
src/SplinesCinterface.o: submodules/Utils/src/Utils/Console.hxx
src/SplinesCinterface.o: submodules/Utils/src/Utils/Malloc.hxx
src/SplinesCinterface.o: submodules/Utils/src/Utils/Numbers.hxx
src/SplinesCinterface.o: submodules/Utils/src/Utils/TicToc.hxx
src/SplinesCinterface.o: submodules/Utils/src/Utils/ThreadPool.hxx
src/SplinesCinterface.o: submodules/GenericContainer/src/GenericContainer.hh
src/SplinesCinterface.o: submodules/GenericContainer/src/GenericContainerConfig.hh
src/SplinesCinterface.o: src/SplinesCinterface.h
src/SplinesUtils.o: src/SplinesUtils.hh src/Splines.hh src/SplinesConfig.hh
src/SplinesUtils.o: submodules/Utils/src/Utils.hh
src/SplinesUtils.o: submodules/Utils/src/Utils/Utils.hxx
src/SplinesUtils.o: submodules/Utils/src/Utils/rang.hpp
src/SplinesUtils.o: submodules/Utils/src/Utils/fmt/printf.h
src/SplinesUtils.o: submodules/Utils/src/Utils/fmt/ostream.h
src/SplinesUtils.o: submodules/Utils/src/Utils/fmt/format.h
src/SplinesUtils.o: submodules/Utils/src/Utils/fmt/core.h
src/SplinesUtils.o: submodules/Utils/src/Utils/fmt/chrono.h
src/SplinesUtils.o: submodules/Utils/src/Utils/fmt/locale.h
src/SplinesUtils.o: submodules/Utils/src/Utils/fmt/ostream.h
src/SplinesUtils.o: submodules/Utils/src/Utils/zstream/izstream.hpp
src/SplinesUtils.o: submodules/Utils/src/Utils/zstream/zstream_common.hpp
src/SplinesUtils.o: submodules/Utils/src/Utils/zstream/izstream_impl.hpp
src/SplinesUtils.o: submodules/Utils/src/Utils/zstream/izstream.hpp
src/SplinesUtils.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
src/SplinesUtils.o: submodules/Utils/src/Utils/zstream/ozstream_impl.hpp
src/SplinesUtils.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
src/SplinesUtils.o: submodules/Utils/src/Utils/Trace.hxx
src/SplinesUtils.o: submodules/Utils/src/Utils/Console.hxx
src/SplinesUtils.o: submodules/Utils/src/Utils/Malloc.hxx
src/SplinesUtils.o: submodules/Utils/src/Utils/Numbers.hxx
src/SplinesUtils.o: submodules/Utils/src/Utils/TicToc.hxx
src/SplinesUtils.o: submodules/Utils/src/Utils/ThreadPool.hxx
src/SplinesUtils.o: submodules/GenericContainer/src/GenericContainer.hh
src/SplinesUtils.o: submodules/GenericContainer/src/GenericContainerConfig.hh
submodules/Utils/src/Console.o: submodules/Utils/src/Utils.hh
submodules/Utils/src/Console.o: submodules/Utils/src/Utils/Utils.hxx
submodules/Utils/src/Console.o: submodules/Utils/src/Utils/rang.hpp
submodules/Utils/src/Console.o: submodules/Utils/src/Utils/fmt/printf.h
submodules/Utils/src/Console.o: submodules/Utils/src/Utils/fmt/ostream.h
submodules/Utils/src/Console.o: submodules/Utils/src/Utils/fmt/format.h
submodules/Utils/src/Console.o: submodules/Utils/src/Utils/fmt/core.h
submodules/Utils/src/Console.o: submodules/Utils/src/Utils/fmt/chrono.h
submodules/Utils/src/Console.o: submodules/Utils/src/Utils/fmt/locale.h
submodules/Utils/src/Console.o: submodules/Utils/src/Utils/fmt/ostream.h
submodules/Utils/src/Console.o: submodules/Utils/src/Utils/zstream/izstream.hpp
submodules/Utils/src/Console.o: submodules/Utils/src/Utils/zstream/zstream_common.hpp
submodules/Utils/src/Console.o: submodules/Utils/src/Utils/zstream/izstream_impl.hpp
submodules/Utils/src/Console.o: submodules/Utils/src/Utils/zstream/izstream.hpp
submodules/Utils/src/Console.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
submodules/Utils/src/Console.o: submodules/Utils/src/Utils/zstream/ozstream_impl.hpp
submodules/Utils/src/Console.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
submodules/Utils/src/Console.o: submodules/Utils/src/Utils/Trace.hxx
submodules/Utils/src/Console.o: submodules/Utils/src/Utils/Console.hxx
submodules/Utils/src/Console.o: submodules/Utils/src/Utils/Malloc.hxx
submodules/Utils/src/Console.o: submodules/Utils/src/Utils/Numbers.hxx
submodules/Utils/src/Console.o: submodules/Utils/src/Utils/TicToc.hxx
submodules/Utils/src/Console.o: submodules/Utils/src/Utils/ThreadPool.hxx
submodules/Utils/src/Malloc.o: submodules/Utils/src/Utils.hh
submodules/Utils/src/Malloc.o: submodules/Utils/src/Utils/Utils.hxx
submodules/Utils/src/Malloc.o: submodules/Utils/src/Utils/rang.hpp
submodules/Utils/src/Malloc.o: submodules/Utils/src/Utils/fmt/printf.h
submodules/Utils/src/Malloc.o: submodules/Utils/src/Utils/fmt/ostream.h
submodules/Utils/src/Malloc.o: submodules/Utils/src/Utils/fmt/format.h
submodules/Utils/src/Malloc.o: submodules/Utils/src/Utils/fmt/core.h
submodules/Utils/src/Malloc.o: submodules/Utils/src/Utils/fmt/chrono.h
submodules/Utils/src/Malloc.o: submodules/Utils/src/Utils/fmt/locale.h
submodules/Utils/src/Malloc.o: submodules/Utils/src/Utils/fmt/ostream.h
submodules/Utils/src/Malloc.o: submodules/Utils/src/Utils/zstream/izstream.hpp
submodules/Utils/src/Malloc.o: submodules/Utils/src/Utils/zstream/zstream_common.hpp
submodules/Utils/src/Malloc.o: submodules/Utils/src/Utils/zstream/izstream_impl.hpp
submodules/Utils/src/Malloc.o: submodules/Utils/src/Utils/zstream/izstream.hpp
submodules/Utils/src/Malloc.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
submodules/Utils/src/Malloc.o: submodules/Utils/src/Utils/zstream/ozstream_impl.hpp
submodules/Utils/src/Malloc.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
submodules/Utils/src/Malloc.o: submodules/Utils/src/Utils/Trace.hxx
submodules/Utils/src/Malloc.o: submodules/Utils/src/Utils/Console.hxx
submodules/Utils/src/Malloc.o: submodules/Utils/src/Utils/Malloc.hxx
submodules/Utils/src/Malloc.o: submodules/Utils/src/Utils/Numbers.hxx
submodules/Utils/src/Malloc.o: submodules/Utils/src/Utils/TicToc.hxx
submodules/Utils/src/Malloc.o: submodules/Utils/src/Utils/ThreadPool.hxx
submodules/Utils/src/Numbers.o: submodules/Utils/src/Utils.hh
submodules/Utils/src/Numbers.o: submodules/Utils/src/Utils/Utils.hxx
submodules/Utils/src/Numbers.o: submodules/Utils/src/Utils/rang.hpp
submodules/Utils/src/Numbers.o: submodules/Utils/src/Utils/fmt/printf.h
submodules/Utils/src/Numbers.o: submodules/Utils/src/Utils/fmt/ostream.h
submodules/Utils/src/Numbers.o: submodules/Utils/src/Utils/fmt/format.h
submodules/Utils/src/Numbers.o: submodules/Utils/src/Utils/fmt/core.h
submodules/Utils/src/Numbers.o: submodules/Utils/src/Utils/fmt/chrono.h
submodules/Utils/src/Numbers.o: submodules/Utils/src/Utils/fmt/locale.h
submodules/Utils/src/Numbers.o: submodules/Utils/src/Utils/fmt/ostream.h
submodules/Utils/src/Numbers.o: submodules/Utils/src/Utils/zstream/izstream.hpp
submodules/Utils/src/Numbers.o: submodules/Utils/src/Utils/zstream/zstream_common.hpp
submodules/Utils/src/Numbers.o: submodules/Utils/src/Utils/zstream/izstream_impl.hpp
submodules/Utils/src/Numbers.o: submodules/Utils/src/Utils/zstream/izstream.hpp
submodules/Utils/src/Numbers.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
submodules/Utils/src/Numbers.o: submodules/Utils/src/Utils/zstream/ozstream_impl.hpp
submodules/Utils/src/Numbers.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
submodules/Utils/src/Numbers.o: submodules/Utils/src/Utils/Trace.hxx
submodules/Utils/src/Numbers.o: submodules/Utils/src/Utils/Console.hxx
submodules/Utils/src/Numbers.o: submodules/Utils/src/Utils/Malloc.hxx
submodules/Utils/src/Numbers.o: submodules/Utils/src/Utils/Numbers.hxx
submodules/Utils/src/Numbers.o: submodules/Utils/src/Utils/TicToc.hxx
submodules/Utils/src/Numbers.o: submodules/Utils/src/Utils/ThreadPool.hxx
submodules/Utils/src/Trace.o: submodules/Utils/src/Utils.hh
submodules/Utils/src/Trace.o: submodules/Utils/src/Utils/Utils.hxx
submodules/Utils/src/Trace.o: submodules/Utils/src/Utils/rang.hpp
submodules/Utils/src/Trace.o: submodules/Utils/src/Utils/fmt/printf.h
submodules/Utils/src/Trace.o: submodules/Utils/src/Utils/fmt/ostream.h
submodules/Utils/src/Trace.o: submodules/Utils/src/Utils/fmt/format.h
submodules/Utils/src/Trace.o: submodules/Utils/src/Utils/fmt/core.h
submodules/Utils/src/Trace.o: submodules/Utils/src/Utils/fmt/chrono.h
submodules/Utils/src/Trace.o: submodules/Utils/src/Utils/fmt/locale.h
submodules/Utils/src/Trace.o: submodules/Utils/src/Utils/fmt/ostream.h
submodules/Utils/src/Trace.o: submodules/Utils/src/Utils/zstream/izstream.hpp
submodules/Utils/src/Trace.o: submodules/Utils/src/Utils/zstream/zstream_common.hpp
submodules/Utils/src/Trace.o: submodules/Utils/src/Utils/zstream/izstream_impl.hpp
submodules/Utils/src/Trace.o: submodules/Utils/src/Utils/zstream/izstream.hpp
submodules/Utils/src/Trace.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
submodules/Utils/src/Trace.o: submodules/Utils/src/Utils/zstream/ozstream_impl.hpp
submodules/Utils/src/Trace.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
submodules/Utils/src/Trace.o: submodules/Utils/src/Utils/Trace.hxx
submodules/Utils/src/Trace.o: submodules/Utils/src/Utils/Console.hxx
submodules/Utils/src/Trace.o: submodules/Utils/src/Utils/Malloc.hxx
submodules/Utils/src/Trace.o: submodules/Utils/src/Utils/Numbers.hxx
submodules/Utils/src/Trace.o: submodules/Utils/src/Utils/TicToc.hxx
submodules/Utils/src/Trace.o: submodules/Utils/src/Utils/ThreadPool.hxx
submodules/Utils/src/Utils.o: submodules/Utils/src/Utils.hh
submodules/Utils/src/Utils.o: submodules/Utils/src/Utils/Utils.hxx
submodules/Utils/src/Utils.o: submodules/Utils/src/Utils/rang.hpp
submodules/Utils/src/Utils.o: submodules/Utils/src/Utils/fmt/printf.h
submodules/Utils/src/Utils.o: submodules/Utils/src/Utils/fmt/ostream.h
submodules/Utils/src/Utils.o: submodules/Utils/src/Utils/fmt/format.h
submodules/Utils/src/Utils.o: submodules/Utils/src/Utils/fmt/core.h
submodules/Utils/src/Utils.o: submodules/Utils/src/Utils/fmt/chrono.h
submodules/Utils/src/Utils.o: submodules/Utils/src/Utils/fmt/locale.h
submodules/Utils/src/Utils.o: submodules/Utils/src/Utils/fmt/ostream.h
submodules/Utils/src/Utils.o: submodules/Utils/src/Utils/zstream/izstream.hpp
submodules/Utils/src/Utils.o: submodules/Utils/src/Utils/zstream/zstream_common.hpp
submodules/Utils/src/Utils.o: submodules/Utils/src/Utils/zstream/izstream_impl.hpp
submodules/Utils/src/Utils.o: submodules/Utils/src/Utils/zstream/izstream.hpp
submodules/Utils/src/Utils.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
submodules/Utils/src/Utils.o: submodules/Utils/src/Utils/zstream/ozstream_impl.hpp
submodules/Utils/src/Utils.o: submodules/Utils/src/Utils/zstream/ozstream.hpp
submodules/Utils/src/Utils.o: submodules/Utils/src/Utils/Trace.hxx
submodules/Utils/src/Utils.o: submodules/Utils/src/Utils/Console.hxx
submodules/Utils/src/Utils.o: submodules/Utils/src/Utils/Malloc.hxx
submodules/Utils/src/Utils.o: submodules/Utils/src/Utils/Numbers.hxx
submodules/Utils/src/Utils.o: submodules/Utils/src/Utils/TicToc.hxx
submodules/Utils/src/Utils.o: submodules/Utils/src/Utils/ThreadPool.hxx
submodules/quarticRootsFlocke/src/PolynomialRoots-1-Quadratic.o: submodules/quarticRootsFlocke/src/PolynomialRoots.hh
submodules/quarticRootsFlocke/src/PolynomialRoots-2-Cubic.o: submodules/quarticRootsFlocke/src/PolynomialRoots.hh
submodules/quarticRootsFlocke/src/PolynomialRoots-3-Quartic.o: submodules/quarticRootsFlocke/src/PolynomialRoots.hh
submodules/quarticRootsFlocke/src/PolynomialRoots-Jenkins-Traub.o: submodules/quarticRootsFlocke/src/PolynomialRoots.hh
submodules/quarticRootsFlocke/src/PolynomialRoots-Utils.o: submodules/quarticRootsFlocke/src/PolynomialRoots.hh
submodules/quarticRootsFlocke/src/PolynomialRoots-Utils.o: submodules/quarticRootsFlocke/src/PolynomialRoots-Utils.hh
submodules/GenericContainer/src/GenericContainer.o: submodules/GenericContainer/src/GenericContainer.hh
submodules/GenericContainer/src/GenericContainer.o: submodules/GenericContainer/src/GenericContainerConfig.hh
submodules/GenericContainer/src/GenericContainerCinterface.o: submodules/GenericContainer/src/GenericContainer.hh
submodules/GenericContainer/src/GenericContainerCinterface.o: submodules/GenericContainer/src/GenericContainerConfig.hh
submodules/GenericContainer/src/GenericContainerCinterface.o: submodules/GenericContainer/src/GenericContainerCinterface.h
submodules/GenericContainer/src/GenericContainerSupport.o: submodules/GenericContainer/src/GenericContainer.hh
submodules/GenericContainer/src/GenericContainerSupport.o: submodules/GenericContainer/src/GenericContainerConfig.hh
submodules/GenericContainer/src/GenericContainerTables.o: submodules/GenericContainer/src/GenericContainer.hh
submodules/GenericContainer/src/GenericContainerTables.o: submodules/GenericContainer/src/GenericContainerConfig.hh
