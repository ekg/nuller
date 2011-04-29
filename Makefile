CXX=g++
BAMTOOLS_ROOT=bamtools
CXXFLAGS=-O3 -I$(BAMTOOLS_ROOT)/include -L./

all: nuller

# builds bamtools static lib, and copies into root
libbamtools.a:
	cd $(BAMTOOLS_ROOT) && mkdir -p build && cd build && cmake .. && make
	cp bamtools/lib/libbamtools.a ./

fastahack/split.o:
	cd fastahack && make

# statically compiles bamaddrg against bamtools static lib
nuller: nuller.cpp libbamtools.a \
	fastahack/Fasta.cpp fastahack/Fasta.h \
	fastahack/split.o \
	vcflib/Variant.cpp vcflib/Variant.h \
	convert.h
	$(CXX) $(CXXFLAGS) nuller.cpp fastahack/split.o fastahack/Fasta.cpp vcflib/Variant.cpp -o nuller -static -lbamtools -lz

clean:
	rm -f nuller libbamtools.a *.o

.PHONY: clean
