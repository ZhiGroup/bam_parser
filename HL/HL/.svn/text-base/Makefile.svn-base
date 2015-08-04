
# current version of the software 
VERSION=0.1

# default installation directory 
INSTALLDIR=/usr/local/bin

# default C++ compiler 
CXX=g++

# default compilation flags are 
# 

SPARSELIB_DIR=../SparseLib++/1.7/
SAMTOOLS_DIR=../samtools-0.1.18
#CFLAGS=-O2 -I../SparseLib++/1.7/include/ -I../SparseLib++/1.7/mv/include/ -I$(SAMTOOLS_DIR) -I./
#CCCFLAGS 	=  -ansi -g -Wall -pedantic -O3 '-DCOMPLEX=std::complex<double>'  
CFLAGS=-I../SparseLib++/1.7/include/ -I../SparseLib++/1.7/mv/include/ -I$(SAMTOOLS_DIR) -I./
#- for debug
CCCFLAGS 	=  -O0 -ansi -ggdb -Wall -pedantic '-DCOMPLEX=std::complex<double>' 
# CCCFLAGS 	=  -O3 -ansi -g -Wall -pedantic '-DCOMPLEX=std::complex<double>'

LDFLAGS	= $(SPARSELIB_DIR)/lib/libsparse.a $(SPARSELIB_DIR)/lib/libspblas.a \
$(SPARSELIB_DIR)/lib/libmv.a -lm -lbam 


# helpful screen listing available options 
help :
	@echo "HL Source Distribution"
	@echo " "
	@echo "This Makefile will compile and install HL 1.0 on your system"
	@echo " "
	@echo "  Type...          To..."
	@echo "  make all         Compile HL"
	@echo "  make clean       Delete temporary files"


bam_reader : bam_reader.o HL.o MyHash.o
	$(CXX) $(CFLAGS) -o $@ bam_reader.o HL.o MyHash.o -L$(SAMTOOLS_DIR) $(LDFLAGS) -lz

bam_reader.o : bam_reader.cpp
	$(CXX) -c $(CFLAGS) $(CCCFLAGS) bam_reader.cpp -o bam_reader.o

bam_parser : bam_parser.o HL.o MyHash.o
	$(CXX) $(CFLAGS) -o $@ bam_parser.o HL.o MyHash.o -L$(SAMTOOLS_DIR) $(LDFLAGS) -static -lz

bam_parser.o : bam_parser.cpp
	$(CXX) -c $(CFLAGS) $(CCCFLAGS) bam_parser.cpp -o bam_parser.o

HL.o : HL.cpp
	$(CXX) -c $(CFLAGS) $(CCCFLAGS) HL.cpp -o HL.o

MyHash.o : MyHash.cpp
	$(CXX) -c $(CFLAGS) $(CCCFLAGS) MyHash.cpp -o MyHash.o

#bam_parser.o : bam_parser.cpp
#	$(CXX) $(CFLAGS) -c $@ bam_parser.o -lz

bgzip:bgzip.o bgzf.o $(KNETFILE_O)
	$(CC) $(CFLAGS) -o $@ bgzf.o bgzip.o $(KNETFILE_O) -lz


EXECUTABLES=bam_parser bam_reader

# make everything 
all : $(EXECUTABLES)

clean :
	rm -f *.o */*.a */*.o $(EXECUTABLES)


.SUFFIXES : .cpp .c .o .X.o $(SUFFIXES)

	
# all: HL sp samtools

# samtools:
# cd ./samtools; make;
	

sp:
	cd ./src; make;
	cd ./spblas; make;
	cd ./mv/src; make ;
	
done:
	@echo "  "
	@echo " Done "
	@echo "  "
