SRCS = \
srcs/SplineAkima.cc \
srcs/SplineAkima2D.cc \
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
srcs/Splines.cc \
srcs/SplinesBivariate.cc \
srcs/SplinesCinterface.cc \
srcs/SplinesUnivariate.cc

OBJS = $(SRCS:.cc=.o)
DEPS = srcs/Splines.hh srcs/SplinesCinterface.h

#CC     = llvm-gcc
#CXX    = llvm-g++
CC     = clang
CXX    = clang++
#CC     = gcc
#CXX    = g++


ifeq ($(shell uname),Linux)
  LIBS   = -static -L./lib -lSplines
  CFLAGS = -I./srcs -I./srcs_utils -Wall -O3 -fPIC -Wno-sign-compare
  AR     = ar rcs
else
  LIBS   = -L./lib -lSplines
  CFLAGS = -I./srcs -I./srcs_utils -Wall -O3 -fPIC -Wno-sign-compare
  AR     = libtool -static -o
endif

all: libSplines.a
	mkdir -p bin
	$(CXX) $(CFLAGS) -o bin/test1 tests/test1.cc $(LIBS)
	$(CXX) $(CFLAGS) -o bin/test2 tests/test2.cc $(LIBS)
	$(CXX) $(CFLAGS) -o bin/test3 tests/test3.cc $(LIBS)
	$(CXX) $(CFLAGS) -o bin/test4 tests/test4.cc $(LIBS)
	$(CXX) $(CFLAGS) -o bin/test5 tests/test5.cc $(LIBS)
	#$(CXX) $(CFLAGS) -o bin/test6 tests/test6.cc $(LIBS)

srcs/%.o: srcs/%.cc $(DEPS)
	$(CXX) $(CFLAGS) -c $< -o $@ 

srcs/%.o: srcs/%.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

libSplines.a: $(OBJS)
	$(AR) lib/libSplines.a $(OBJS) 

libSplines.dylib: $(OBJS)
	$(CXX) -shared -o lib/libSplines.dylib $(OBJS) 

libSplines.so: $(OBJS)
	$(CXX) -shared -o lib/libSplines.so $(OBJS) 

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
	rm -f lib/libSplines.a lib/libSplines.dylib lib/libSplines.so srcs/*.o
	rm -f tests/test1
	rm -rf bin
	