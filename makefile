# get the type of OS currently running
OS=$(shell uname)
PWD=$(shell pwd)

LIB_SPLINE = libSplines.a
LIB_GC     = libGenericContainer.a

CC   = gcc
CXX  = g++
INC  = -I./srcs -I./include -I./GC/srcs
LIBS = -L./lib -lSplines -lGenericContainer
DEFS =

# check if the OS string contains 'Linux'
ifneq (,$(findstring Linux, $(OS)))
  #LIB_SPLINE = libSplines.so
  #LIB_GC     = libGenericContainer.so
  LIBS   = -static -L./lib -lSplines -lGenericContainer
  CFLAGS = -Wall -O3 -fPIC -Wno-sign-compare
  AR     = ar rcs
endif

# check if the OS string contains 'Darwin'
ifneq (,$(findstring Darwin, $(OS)))
  CC     = clang
  CXX    = clang++
  #LIB_SPLINE = libSplines.dylib
  #LIB_GC     = libGenericContainer.dylib
  LIBS   = -L./lib -lSplines -lGenericContainer
  CFLAGS = -Wall -O3 -fPIC -Wno-sign-compare
  AR     = libtool -static -o
endif

SRCS = \
srcs/SplineAkima.cc \
srcs/SplineAkima2D.cc \
srcs/SplineBSpline.cc \
srcs/SplineBessel.cc \
srcs/SplineBiCubic.cc \
srcs/SplineBiQuintic.cc \
srcs/SplineBilinear.cc \
srcs/SplineConstant.cc \
srcs/SplineCubic.cc \
srcs/SplineCubicBase.cc \
srcs/SplineHermite.cc \
srcs/SplineLinear.cc \
srcs/SplinePchip.cc \
srcs/SplineQuintic.cc \
srcs/SplineQuinticBase.cc \
srcs/SplineSet.cc \
srcs/SplineSetGC.cc \
srcs/Splines.cc \
srcs/SplinesBivariate.cc \
srcs/SplinesCinterface.cc \
srcs/SplinesUnivariate.cc

OBJS  = $(SRCS:.cc=.o)
DEPS  = srcs/Splines.hh srcs/SplinesCinterface.h
MKDIR = mkdir -p

# prefix for installation, use make PREFIX=/new/prefix install
# to override
PREFIX    = /usr/local
FRAMEWORK = Splines

all: lib/$(LIB_SPLINE) lib/$(LIB_GC)
	mkdir -p bin
	$(CXX) $(INC) $(CFLAGS) -o bin/test1 tests/test1.cc $(LIBS)
	$(CXX) $(INC) $(CFLAGS) -o bin/test2 tests/test2.cc $(LIBS)
	$(CXX) $(INC) $(CFLAGS) -o bin/test3 tests/test3.cc $(LIBS)
	$(CXX) $(INC) $(CFLAGS) -o bin/test4 tests/test4.cc $(LIBS)
	$(CXX) $(INC) $(CFLAGS) -o bin/test5 tests/test5.cc $(LIBS)
	#$(CXX) $(CFLAGS) -o bin/test6 tests/test6.cc $(LIBS)

srcs/%.o: srcs/%.cc $(DEPS)
	$(CXX) $(INC) $(CFLAGS) $(DEFS) -c $< -o $@ 

srcs/%.o: srcs/%.c $(DEPS)
	$(CC) $(INC) $(DEFS) -c -o $@ $< $(CFLAGS)

lib/libSplines.a: $(OBJS)
	$(AR) lib/libSplines.a $(OBJS) 

lib/libSplines.dylib: $(OBJS)
	$(CXX) -shared -o lib/libSplines.dylib $(OBJS) 

lib/libSplines.so: $(OBJS)
	$(CXX) -shared -o lib/libSplines.so $(OBJS) 

lib/$(LIB_GC):
	$(MKDIR) include ; cd GC ; make ; make PREFIX=$(PWD) install 


install: lib/$(LIB_SPLINE) lib/$(LIB_GC)
	cp srcs/Splines.hh          $(PREFIX)/include
	cp srcs/SplinesCinterface.h $(PREFIX)/include
	cp lib/$(LIB_SPLINE)        $(PREFIX)/lib

install_as_framework: lib/$(LIB_GC)
	$(MKDIR) $(PREFIX)/include/$(FRAMEWORK)
	cp srcs/Splines.hh          $(PREFIX)/include/$(FRAMEWORK)
	cp srcs/SplinesCinterface.h $(PREFIX)/include/$(FRAMEWORK)
	cp lib/$(LIB_SPLINE)        $(PREFIX)/lib

run:
	./bin/test1
	./bin/test2
	./bin/test3
	./bin/test4
	./bin/test5
	#./bin/test6

doc:
	doxygen
	
clean:
	rm -f lib/libSplines.* lib/libGenericContainer.* srcs/*.o

	rm -rf bin
	