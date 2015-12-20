.PHONY: readSim
all: readSim


# I am a comment, and I want to say that the variable CC will be
# the compiler to use.
CC=g++
# Hey!, I am comment number 2. I want to say that CFLAGS will be the
# options I'll pass to the compiler.
CFLAGS=-c -Wall

htslib/libhts.so:
	git clone https://github.com/samtools/htslib.git
	cd htslib;make

readSim: htslib/libhts.so 
	 $(CC)  -g -I src/ -Wall -c -o src/empdist.o src/empdist.cpp
	 $(CC)  -g -std=c++0x -I htslib/htslib -L htslib -I /r1/people/homa_amini/Repository/SimulateGenome/boost_1_57_0 -Wl,-rpath=/r1/people/homa_amini/Repository/SimulateGenome/boost_1_57_0/stage/lib/ -Wl,-rpath=/r1/people/homa_amini/Repository/SimulateGenome/htslib src/readSim.cpp -o readSim src/empdist.o -lhts -lpopt -pthread -lboost_random 
	# $(CC)  -g -std=c++0x -I /r1/people/homa_amini/Repository/SimulateGenome/boost_1_57_0 -Wl,-rpath=/r1/people/homa_amini/Repository/SimulateGenome/boost_1_57_0/stage/lib/ -Wl,-rpath=/r1/people/homa_amini/Repository/SimulateGenome/htslib src/readSim.cpp -o readSim src/empdist.o -lhts -lpopt -pthread -lboost_random 

readSimDip: htslib/libhts.so 
	 $(CC)  -g -I src/ -Wall -c -o src/empdist.o src/empdist.cpp
	 $(CC)  -g -std=c++0x -I htslib/htslib -L htslib -I /r1/people/homa_amini/Repository/SimulateGenome/boost_1_57_0 -Wl,-rpath=/r1/people/homa_amini/Repository/SimulateGenome/boost_1_57_0/stage/lib/ -Wl,-rpath=/r1/people/homa_amini/Repository/SimulateGenome/htslib src/readSimDiploid.cpp -o readSimDiploid src/empdist.o -lhts -lpopt -pthread -lboost_random 
	# $(CC)  -g -std=c++0x -I /r1/people/homa_amini/Repository/SimulateGenome/boost_1_57_0 -Wl,-rpath=/r1/people/homa_amini/Repository/SimulateGenome/boost_1_57_0/stage/lib/ -Wl,-rpath=/r1/people/homa_amini/Repository/SimulateGenome/htslib src/readSim.cpp -o readSim src/empdist.o -lhts -lpopt -pthread -lboost_random 

install:
	install readSim /usr/local/bin/readSim

clean:
	rm -rf src/readSim.o readSim -rf htslib
