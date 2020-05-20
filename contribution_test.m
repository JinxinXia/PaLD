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
d = rand(1,5,5);
d
C = get_ctrib_mat(d,1);
C









