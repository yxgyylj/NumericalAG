%================ Simple check for residues of fineGrid_bdry =============

N = 160;
load('Data/parameters.mat');
load(sprintf('Data/solution_brdy_mu_2_%d.mat',N));

% print residues
res = zeros(1,size(solutions,2));
for j = 1:size(solutions,2)
    foo = GS_RHS_brdy(solutions(:,j),p);
    res(j) = norm(foo);
end
TOL = 1e-6;
good_ind = find(res<TOL);
bad_ind = find(res>=TOL);