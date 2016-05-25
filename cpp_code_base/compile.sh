#!/usr/bin/sh

rm build/*.o
rm model_sim

g++ -Wall -c -o build/solver.o source/solver.cpp
g++ -Wall -c -o build/main.o source/main.cpp
g++ -Wall -c -o build/integrator.o source/integrator.cpp
g++ -Wall -c -o build/contactdetector.o source/contactdetector.cpp
g++ -Wall -c -o build/state.o source/state.cpp
g++ -Wall -c -o build/parameters.o source/parameters.cpp
g++ -Wall -o model_sim build/parameters.o build/main.o build/integrator.o build/contactdetector.o build/solver.o build/state.o