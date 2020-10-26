
clear;
rng(4);

% create a random distance matrix that is symmetric with diagonal elements
% equal to zeros
d = rand(1000);
D = (d+d')/2;
D = D - diag(diag(D));

profile on

disp('their method')

tic
[C1,F] = orig_contribute(D,1);
toc


% b is different based on cache size, the machine cache size 
% in my computer is 9MB, which is around the size that could 
% store 1000 x 1000 matrix

% in the raw block method, the b should equal to sqrt(M)/3
% where M is the cache size
disp('block method')
tic
[C2,U] = matmul_block(D,1,500);
toc 

disp('optimal block method')
tic
[C3,U1] = optimal_matmul_block(D,1,500);
toc 


norm(C1-C3)





