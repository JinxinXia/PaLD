CC = gcc
CFLAGS = -Wall -O3
LDFLAGS = -lm

objects = PaLD_test.o kernels.o utils.o
targets = PaLD_test

.PHONY : default
default: $(targets)

.PHONY : all
all: clean $(targets)

PaLD_test : $(objects)
	$(CC) $^ -o $@ $(LDFLAGS) 

%.o : %.c %.h
	$(CC) -c $(CFLAGS) $< -o $@

.PHONY : clean
clean :
	rm -rf $(objects) $(targets) *~


