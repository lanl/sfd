#The C compiler to use
CPP = /usr/bin/g++
# OPEN MP aware compiler on OS X
#CPP = /usr/local/bin/g++-6

#Set your path to the executable
MYPATH = $(PWD)

#Code version
PROG = Moment
VERSION = 5

# OPENMP
#CFLAGS = -O3 -fopenmp -DOMP
# No OPENMP
CFLAGS = -O3

SRC_FILES = *.cpp
OBJ = *.o
MYLIB = $(PROG).$(VERSION).a

default: main

# Compile the code
main:lib
	$(CPP) -o ./$(PROG).$(VERSION) ./Main.cpp -lm $(CFLAGS) $(MYLIB)

lib:	$(PROG).$(VERSION).a

$(PROG).$(VERSION).a:	$(OBJ)
	ar ruv $(PROG).$(VERSION).a $(OBJ)
	ranlib $(PROG).$(VERSION).a

$(OBJ):
	$(CPP) $(CFLAGS) -c $(SRC_FILES)

# Clean up some stuff
clean:
	rm -f *.o

clean_all:
	rm -f *.o ./$(PROG).$(VERSION) *.a ./Moment

wrapper:
	echo "#!/bin/csh\n$(MYPATH)/$(PROG).$(VERSION) | tee -a log.out" > ./Moment
	chmod +x ./Moment


test:
	/usr/bin/time ./$(PROG).$(VERSION)












