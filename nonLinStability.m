%=================== Iterative stability of 1d Problem ====================

clear; close all; clc;
addpath(genpath('src'));

%% Parameters
load('Data/solution_160.mat');
load('Data/parameters.mat');
nGrid = floor(size(solution,1)/2);
xx = linspace(p.xl,p.xr,nGrid)';
nSol = size(solution,2);     % # of subdomains
tspan = 0:.02:10;
tolerence = 1e-2;
ind = [];

%% Start iteration
tic;
fprintf('There are in total %d solutions. \n', nSol);
for k = 1:nSol
    fprintf('Checking for the %d-th solution. \n', k);
    noise = .01*mean(abs(solution(:,k)))*(rand(size(solution,1),1)-.5);
    guess = solution(:,k) + noise;
    [tt, sol_temp] = ode15s(@(t,x)(GS_RHS_t(t,x,p)),tspan,guess);
    
    %% Display solutions
    subplot(3,4,k);
    plot(xx,solution(1:nGrid,k))
    hold on;
    plot(xx,sol_temp(end,1:nGrid)')
    plot(xx,guess(1:nGrid)')
    xlabel('x');
    ylabel('A');
    legend('SSS','Iterated','Guess')
    title(sprintf('The %d-th solution', k));
    
    if norm(sol_temp(end,:)' - solution(:,k)) < tolerence
        ind = [ind, k];
    end
end

savefig( sprintf('Imgs/Nonlinear_test_%d.fig', nGrid-1));
end_time = toc;
fprintf('Elapsed time = %f seconds \n', end_time);

% %% Plot final outcome
% figure(1);
% set(gca,'fontsize',20);
% subplot(2,1,1)
% plot(xx,solution(1:nGrid,ind))
% xlabel('x');
% ylabel('A');
% title(sprintf('Solutions in N=%d subdomains', 160));
% subplot(2,1,2)
% plot(xx,solution(nGrid+1:2*nGrid,ind))
% xlabel('x');
% ylabel('S');
% title(sprintf('Number of stable solutions = %d', length(ind)));