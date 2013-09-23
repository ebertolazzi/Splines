SRCS = \
srcs/SplineAkima.cc \
srcs/SplineBessel.cc \
srcs/SplineCubic.cc \
srcs/SplineCubicBase.cc \
srcs/SplineHermite.cc \
srcs/SplinePchip.cc \
srcs/SplineQuintic.cc \
srcs/SplineQuinticBase.cc \
srcs/Splines.cc

OBJS = $(SRCS:.cc=.o)
DEPS = srcs/Splines.hh

#CC     = llvm-gcc
#CXX    = llvm-g++
#CC     = clang
#CXX    = clang++
CC     = gcc
CXX    = g++

CFLAGS =  -I./srcs -Wall -O3
LIBS   = -L./lib -lSplines

#AR     = ar rcs
AR     = libtool -static -o 

all: libSplines.a
	$(CXX) $(CFLAGS) -o bin/test1 tests/test1.cc $(LIBS)

srcs/%.o: srcs/%.cc $(DEPS)
	$(CXX) $(CFLAGS) -c $< -o $@ 

srcs/%.o: srcs/%.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

libSplines.a: $(OBJS)
	$(AR) lib/libSplines.a $(OBJS) 

run:
	cd bin ; ./test1

doc:
	doxygen
	
clean:
	rm -f lib/libSplines.a srcs/*.o
	rm -f tests/test1
	