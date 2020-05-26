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
d = rand(500,500);
D = (d+d')/2;
D = D - diag(diag(D));

disp('our method')

tic
C1 = get_ctrib_mat(D,1);
toc

disp('their method')
tic
C = get_contribution_matrix(D,1);
toc
fprintf('norm(C1-C): %f \n',norm(C-C1))
