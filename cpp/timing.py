import os

start = 500
end   = 5000
for cache_size in [8, 16, 32, 64, 128]:
    print('cache_size is ' + str(cache_size))
    for matrix_size in range(start, end, int((end-start-1)/49)):
        os.system('./PaLD_test ' + str(matrix_size) + ' ' + str(cache_size) + 
                ' >> exp_cache_size_' + str(cache_size)+'.txt')
