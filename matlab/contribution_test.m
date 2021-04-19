clear;
rng(6);

% create a random distance matrix that is symmetric with diagonal elements
% equal to zeros
n = 5; b1 = 100; b2 = 2;
d = rand(n);
%d = [0 1 2 3; 1 0 4 5; 2 4 0 6; 3 5 6 0];
D = (d+d')/2;
D = D - diag(diag(D));

%profile on

fprintf('their method, running on dimension %d\n', n)

tic
[C1,F] = pald_orig(D,1);
toc


% b is different based on cache size, the machine cache size 
% in my computer is 9MB, which is around the size that could 
% store 1000 x 1000 matrix

% in the raw block method, the b should equal to sqrt(M)/3
% where M is the cache size

%{
fprintf('block method, same dimension, using block size %d\n', b1)
tic
[C2,U] = pald_block(D,1,b1);
toc 
fprintf('error in blocked method is %g\n', norm(C1-C2))

fprintf('optimal block method, same dimension, using block size %d\n', b2)
tic
[C3,U1] = pald_opt(D,1,b2);
toc 
fprintf('error in optimal blocked method is %g\n', norm(C1-C3))
%}



fprintf('triplet method, same dimension\n')
tic
[C4,U4] = pald_triplet(D);
toc 
fprintf('error in triplet method is %g\n', norm(C1-C4))

fprintf('triplet block method, same dimension\n')
tic
[C5,U5] = pald_triplet_block(D,b2);
toc 
fprintf('error in triplet method is %g\n', norm(C1-C5))

