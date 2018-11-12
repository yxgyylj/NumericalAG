%==================== Linear stability of 1d Problem ======================

clear; close all; clc;
addpath(genpath('src'));

%% Parameters
load('Data/solution_160.mat');
load('Data/parameters.mat');
nGrid = floor(size(solution,1)/2);
xx = linspace(p.xl,p.xr,nGrid)';
nSol = size(solution,2);     % # of subdomains
ind = [];

%% Main loop -- find the max eigenvalue of the linearized problem
for k = 1:nSol
    J = GS_Jacobian(nGrid,solution(:,k),p);
    if max(eig(J)) < 0
        ind = [ind k];
    end
end

%% Plot final outcome
figure(1);
subplot(2,1,1)
plot(xx,solution(1:nGrid,ind))
xlabel('x');
ylabel('A');
title(sprintf('Solutions in N=%d subdomains', nGrid-1));
set(gca,'fontsize',20);
subplot(2,1,2)
plot(xx,solution(nGrid+1:2*nGrid,ind))
xlabel('x');
ylabel('S');
title(sprintf('Number of stable solutions = %d', length(ind)));
set(gca,'fontsize',20);
saveas(gcf, sprintf('Imgs/Linear_stability_%d.png', nGrid-1));