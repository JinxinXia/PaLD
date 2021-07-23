CC =icc
CFLAGS = -Wall -O3 -qopenmp -xHost -g -pg -DBLOCKSIZE=256 -DBLOCKSIZEPAR=256 -qopt-zmm-usage=high
LDFLAGS = -lm -qopenmp

objects = original.o opt.o utils.o opt_partial.o
targets = PaLD_test PaLD_par_test PaLD_test_avx

.PHONY : default
default: $(targets)

.PHONY : all
all: clean $(targets)

PaLD_test : PaLD_test.o $(objects)
	$(CC) $^ -o $@ $(LDFLAGS) 

PaLD_par_test : PaLD_par_test.o $(objects)
	$(CC) $^ -o $@ $(LDFLAGS) 

PaLD_test_avx : PaLD_test_avx.o $(objects)
	$(CC) $^ -o $@ $(LDFLAGS) 

%.o : %.c %.h
	$(CC) -c $(CFLAGS) $< -o $@

.PHONY : clean
clean :
	rm -rf $(objects) $(targets) *~


