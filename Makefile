CC=c++ 
LDFLAGS =
CFLAGS = -c -fPIC -DDEBUG -g -std=c++11
LIBS=-lntl -lgmp
LIB_DIR = -L/usr/local/lib
INCLUDES = -I/usr/local/include/NTL		
COMMON_SOURCES = ffspol.cpp ffstools.cpp poly_q.cpp ij_vector.cpp q_lattice.cpp ffsnorm.cpp bucket.cpp factorBaseMisc.cpp smoothnessTest.cpp latticeSieve.cpp 
FACTORBASE_SOURCES = factor_base.cpp 
POLYSELECT_SOURCES = polynomial_selection.cpp 
SIEVE_SOURCES= sieve.cpp
COMMON_OBJECTS=$(COMMON_SOURCES:.cpp=.o)
FACTORBASE_OBJECTS=$(FACTORBASE_SOURCES:.cpp=.o)
POLYSELECT_OBJECTS=$(POLYSELECT_SOURCES:.cpp=.o)
SIEVE_OBJECTS=$(SIEVE_SOURCES:.cpp=.o)	
FACTOR_EXECUTABLE = factorBase
POLY_EXECUTABLE = polySelect
SIEVE_EXECUTABLE = sieve

.PHONY: all factorside polyselectside sieveside

all: factorside polyselectside sieveside

factorside: $(FACTOR_EXECUTABLE)

polyselectside: $(POLY_EXECUTABLE)
	
sieveside: $(SIEVE_EXECUTABLE)	

$(FACTOR_EXECUTABLE):$(COMMON_OBJECTS) $(FACTORBASE_OBJECTS)
	$(CC) $(LDFLAGS)  $^ -o $@ $(LIB_DIR) $(LIBS) $(INCLUDES) 

$(POLY_EXECUTABLE):$(COMMON_OBJECTS) $(POLYSELECT_OBJECTS)
	 $(CC) $(LDFLAGS)  $^ -o $@ $(LIB_DIR) $(LIBS) $(INCLUDES) 

$(SIEVE_EXECUTABLE):$(COMMON_OBJECTS) $(SIEVE_OBJECTS)
	 $(CC) $(LDFLAGS)  $^ -o $@ $(LIB_DIR) $(LIBS) $(INCLUDES)	

.cpp.o:
	$(CC) $(CFLAGS) $(LIB_DIR) $(LIBS) $(INCLUDES)  $< -o $@ 

clean:
	rm -rf $(FACTOR_EXECUTABLE) *.o
	rm -rf $(POLY_EXECUTABLE) *.o
	rm -rf $(SIEVE_EXECUTABLE) *.o 