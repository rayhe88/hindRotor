main: main.o Rvector.o Atom.o Molecule.o Inertia.o Matrix.o HPot.o FourierHamiltonian.o thermo.o runCommands.o
	g++ -g -Wall -o hRotor.x main.o Rvector.o Atom.o Molecule.o Inertia.o Matrix.o HPot.o FourierHamiltonian.o thermo.o runCommands.o -llapack -llapack

main.o: main.cpp Rvector.cpp Rvector.h
	g++ -g -Wall -c main.cpp

Rvector.o: Rvector.cpp Rvector.h
	g++ -g -Wall -c Rvector.cpp

Atom.o: Atom.cpp Atom.h
	g++ -g -Wall -c Atom.cpp

Molecule.o: Molecule.cpp Molecule.h screen.h
	g++ -g -Wall -c Molecule.cpp

Inertia.o: Inertia.cpp Inertia.h
	g++ -g -Wall -c Inertia.cpp

Matrix.o: Matrix.cpp Matrix.h
	g++ -g -Wall -c Matrix.cpp

FourierHamiltonian.o: FourierHamiltonian.cpp FourierHamiltonian.h screen.h
	g++ -g -Wall -c FourierHamiltonian.cpp

HPot.o: HPot.cpp HPot.h
	g++ -g -Wall -c HPot.cpp

thermo.o: thermo.cpp thermo.h
	g++ -g -Wall -c thermo.cpp

runCommands.o: runCommands.cpp runCommands.h screen.h
	g++ -g -Wall -c runCommands.cpp
    
.PHONY: clean
clean:
	rm *.o hRotor.x
