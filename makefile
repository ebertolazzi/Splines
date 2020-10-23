# get the type of OS currently running
OS=$(shell uname)
PWD=$(shell pwd)

LIB_SPLINE = libSplines.a
LIB_GC     = libGenericContainer.a

CC   = gcc
CXX  = g++ -std=c++11
INC  = -Isrc -Iinclude -IGC/lib/include -Isubmodules/Utils/src -Isubmodules/quarticRootsFlocke/src
LIBS = -Llib -LGC/lib -lSplines -lGenericContainer
DEFS =
MAKE = make

# check if the OS string contains 'Linux'
ifneq (,$(findstring Linux, $(OS)))
  WARN = -Wall
  # activate C++11 for g++ >= 4.9
  VERSION  = $(shell $(CC) -dumpversion)
  CC      += $(WARN)
  CXX     += $(WARN)
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
        $(shell echo submodules/Utils/src/fmt/*.cc) \
        $(shell echo submodules/quarticRootsFlocke/src/*.cc)
OBJS  = $(SRCS:.cc=.o)
DEPS  = $(shell echo src/*.hh) $(shell echo submodules/Utils/src/*.h*)
MKDIR = mkdir -p

# prefix for installation, use make PREFIX=/new/prefix install
# to override
PREFIX    = /usr/local
FRAMEWORK = Splines

all: gc lib tests

gc: lib/$(LIB_GC)

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

travis: gc lib tests run

include_local: gc
	@rm -rf lib/include
	$(MKDIR) lib
	$(MKDIR) lib/include
	@cp -f src/*.h* lib/include

src/%.o: src/%.cc $(DEPS) gc
	$(CXX) $(INC) $(CXXFLAGS) $(DEFS) -c $< -o $@

submodules/Utils/src/%.o: submodules/Utils/src/%.cc $(DEPS) gc
	$(CXX) $(INC) $(CXXFLAGS) $(DEFS) -c $< -o $@

submodules/Utils/src/fmt/%.o: submodules/Utils/src/fmt/%.cc $(DEPS) gc
	$(CXX) $(INC) $(CXXFLAGS) $(DEFS) -c $< -o $@

src/%.o: src/%.c $(DEPS) gc
	$(CC) $(INC) $(CFLAGS) $(DEFS) -c -o $@ $<

lib/libSplines.a: $(OBJS) include_local
	$(AR) lib/libSplines.a $(OBJS)

lib/libSplines.dylib: $(OBJS) include_local
	$(CXX) -shared -o lib/libSplines.dylib $(OBJS)

lib/libSplines.so: $(OBJS) include_local
	$(CXX) -shared -o lib/libSplines.so $(OBJS)

lib/$(LIB_GC):
	rm -rf GC; git clone -b develop --depth 1 https://github.com/ebertolazzi/GenericContainer.git GC
ifneq (,$(findstring Linux, $(OS)))
	cd GC; ruby gcc_workaround.rb; cd ..
endif
	$(MKDIR) include; cd GC; ${MAKE} CXXFLAGS="$(CXXFLAGS)" CC="$(CC)" CXX="$(CXX)" lib;

install_local: lib
	$(MKDIR) ./lib/include
	cp GC/lib/include/*        ./lib/include
	cp src/Splines.hh          ./lib/include
	cp src/SplinesCinterface.h ./lib/include

install: lib
	cp src/Splines.hh          $(PREFIX)/include
	cp src/SplinesCinterface.h $(PREFIX)/include
	cp lib/$(LIB_SPLINE)       $(PREFIX)/lib

install_as_framework: lib
	$(MKDIR) $(PREFIX)/include/$(FRAMEWORK)
	cp src/Splines.hh          $(PREFIX)/include/$(FRAMEWORK)
	cp src/SplinesCinterface.h $(PREFIX)/include/$(FRAMEWORK)
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
	doxygen

clean:
	rm -rf lib/libSplines.* lib/libGenericContainer.* lib/include src/*.o
	rm -rf bin
