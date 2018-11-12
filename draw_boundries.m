%======================= Draw parameter boundaries ========================
%   A bifurcation diagram can be used to illustrate how the number of 
%   solutions can change as parameters vary. Numerically, drawing a
%   bifurcation diagram can be proceed by calculating the corresponding
%   perturbation problem.
%
%       Step 1: derive the corresponding small perturbation problem,
%       collecting first order terms;
%       Step 2: solving the first order problem along with the original
%       problem;
%       Step 3: draw boundaries
%
%   by Xige Yang
%=========================================================================

clear; close all; clc;
addpath(genpath('src'));

%% Extract solutions
mu_start = 1;
mu_end = 0;
load(sprintf('Data/track_rho_mu_%d_%d_160t', mu_start,mu_end = 0),'solutions','finer_grid');

%% Plot final outcome
figure(1);
subplot(2,2,1)
plot(finer_grid,solutions(1:nGrid,:))
xlabel('x');
ylabel('A');
title(sprintf('Solutions in N=%d subdomains', p.N*2^j));
set(gca,'fontsize',16);

subplot(2,2,3)
plot(finer_grid,solutions(nGrid+1:2*nGrid,:))
xlabel('x');
ylabel('S');
title(sprintf('Number of solutions = %d', nSol_finer));
set(gca,'fontsize',16);

subplot(2,2,2)
plot(finer_grid,solutions(2*nGrid_old+1:3*nGrid_old,:))
xlabel('x');
ylabel('A1');
title('First order terms');
set(gca,'fontsize',16);

subplot(2,2,4)
plot(finer_grid,solutions(3*nGrid_old+1:4*nGrid_old,:))
xlabel('x');
ylabel('S1');
set(gca,'fontsize',16);
savefig(sprintf('Imgs/bdry_mu_%d_N_%d.fig', mu_start, p.N*2^j));

figure(2);
subplot(2,1,1)
plot(1:nSol_finer,solutions(end,:)','b*');
title('rhos');
xlabel('solution id');
ylabel(['rho']);
set(gca,'fontsize',16);

subplot(2,1,2)
plot(1:nSol_finer,res,'r*');
title('Residual plot');
xlabel('solution id');
ylabel(['norm of Newton\prime', 's method']);
set(gca,'fontsize',16);