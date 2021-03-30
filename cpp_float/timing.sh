#!/bin/bash

for matrix_size in 3000 4000 5000
do
    for cache_size in 16 32 64 128 256 512
    do
	./PaLD_test $matrix_size $cache_size >> icc_f_$cache_size.txt
    done
done
