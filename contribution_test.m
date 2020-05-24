clear;
%T = readtable('distance.csv');


%cell_dist = T{:,:};
%b = cell_dist(2:end,2:end);


%class(b)
%size(b)

%B = [b{:,:}];

%class(B)
%C = get_contribution_matrix(b,1);


clear;
rng(12);
d = rand(5,5);
D = (d+d')/2;
D



tic
C1 = get_ctrib_mat(d,1);
toc

tic
C = get_contribution_matrix(d,1);
toc
