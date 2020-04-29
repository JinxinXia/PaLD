clear;
T = readtable('distance.csv');


cell_dist = T{:,:};
b = cell_dist(2:end,2:end);


class(b)
size(b)

%C = get_contribution_matrix(b,1);