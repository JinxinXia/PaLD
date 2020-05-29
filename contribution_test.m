
clear;
rng(122);

% create a random distance matrix that is symmetric with diagonal elements
% equal to zeros
d = rand(50,50);
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


disp('their method')
tic
C2 = get_cmat(D,1);
toc
% check distance between the output of two methods
fprintf('norm(C-C1): %f \n',norm(C-C1))
fprintf('norm(C1-C2): %f \n',norm(C1-C2))
