CXX = mpic++
CC = mpic++
LDFLAGS = -lstdc++
CPPFLAGS = -Ofast
CXXFLAGS = -std=c++17 -fpermissive
LDLIBS = -lboost_timer -lboost_program_options -lm -lgsl -lgslcblas -lga
OBJ0 = Polar.o PolarGenome.o MPIparaMain.o # test.o
# OBJ1 = # accumulate.o
EXE0 = MPIparaMain # test
# EXE1 = # accumulate
OBJECTS = ${OBJ0} #/usr/lib/libga.a # ${OBJ1}
EXECUTABLE = ${EXE0} # ${EXE1}

all: ${EXECUTABLE}
	./GenomeRun --help
${EXECUTABLE}: ${OBJECTS} 
# ${EXE1}: ${OBJ1}
clean:
	rm ${OBJ0} ${EXE0}
