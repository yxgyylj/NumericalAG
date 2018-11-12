%=============== Compare coarse-grid with finer-grid solutions ============
%
%   by Xige Yang
%==========================================================================

clear; close all; clc;
%% Extract solutions
load('Data/parameters_sep.mat');
p.N = 5;
freal='bertiniFile/Sep/input_5/real_finite_solutions';
[init_sol, nSol] = read_real_sol(freal,p.N);

%% Plot original solutions
subplot(3,2,1)
plot(x,init_sol(1:length(x),:));
xlabel('x');
ylabel('A');
ylim([0,.4]);
title('Original A at N = 5');
set(gca,'fontsize',16);

%% Plot refined solutions
for j = 1:5
    N_cur = p.N*2^j;
    load(sprintf('Data/pure_refine_sep_N_5/solutions_%d.mat',N_cur));
    x = finer_grid;
    subplot(3,2,j+1)
    plot(x,solutions(1:length(x),:));
    xlabel('x');
    ylabel('A');
    ylim([0,.4]);
    title(sprintf('N = %d',N_cur));
    set(gca,'fontsize',16);
end

savefig('Imgs/coarse_vs_finer.fig')