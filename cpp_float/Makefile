CC = icc
CFLAGS = -Ofast -Wall -qopenmp -xHost -g -DEBUG
LDFLAGS = -lm -qopenmp -g -DEBUG

objects = kernels.o utils.o
targets = PaLD_test PaLD_par_test

.PHONY : default
default: $(targets)

.PHONY : all
all: clean $(targets)

PaLD_test : PaLD_test.o $(objects)
	$(CC) $^ -o $@ $(LDFLAGS) 

PaLD_par_test : PaLD_par_test.o $(objects)
	$(CC) $^ -o $@ $(LDFLAGS) 

%.o : %.c %.h
	$(CC) -c $(CFLAGS) $< -o $@

.PHONY : clean
clean :
	rm -rf $(objects) $(targets) *~ *.o


