CC = g++
CCFLAG = -std=c++17
INC = $(PWD)

all : testlive.exe

test : test.exe

testlive.exe : makefile testlive.cpp Node.hpp ParaBox.hpp MersenneTwister.h LatticeSite.hpp SquareLattice.hpp Diagram.hpp Worm.hpp WLDiagram.hpp DipoleWorm.hpp LiveDipoleWorm.hpp Check.hpp ReadInput.h NUpdate.hpp Measure.hpp
	$(CC) $(CCFLAG) testlive.cpp -o testlive.exe -I$(INC)

test.exe : makefile test.cpp Node.hpp ParaBox.hpp MersenneTwister.h LatticeSite.hpp SquareLattice.hpp Diagram.hpp Worm.hpp WLDiagram.hpp DipoleWorm.hpp ReadInput.h
	$(CC) $(CCFLAG) test.cpp -o test.exe -I$(INC)

clean :
	rm *.exe *.diag *.wld *~ *.info *.backup
