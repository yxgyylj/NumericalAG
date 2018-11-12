%=============== Compare coarse-grid with finer-grid solutions ============
%
%   by Xige Yang
%==========================================================================

clear; close all; clc;
%% Extract solutions
load('Data/parameters.mat');
freal='bertiniFile/real_finite_solutions';
[init_sol, nSol] = read_real_sol(freal,p.N);

%% Plot original solutions
subplot(2,2,1)
plot(x,init_sol(1:length(x),:));
xlabel('x');
ylabel('A');
title('N = 5');
set(gca,'fontsize',16);
subplot(2,2,3)
plot(x,init_sol(length(x)+1:end,:));
xlabel('x');
ylabel('S');
title(sprintf('# of solution = %d', size(init_sol,2)));
set(gca,'fontsize',16);

%% Plot refined solutions
load('Data/solution_160.mat');
x = finer_grid;
subplot(2,2,2)
plot(x,solution(1:length(x),:));
xlabel('x');
ylabel('A');
title('N = 160');
set(gca,'fontsize',16);
subplot(2,2,4)
plot(x,solution(length(x)+1:end,:));
xlabel('x');
ylabel('S');
title(sprintf('# of solution = %d', size(solution,2)));
set(gca,'fontsize',16);

saveas(gcf,'Imgs/coarse_vs_finer.png')