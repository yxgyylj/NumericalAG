%======================= Test file for track_rho.m ========================
%
%   Step 1: load data (which is got by fineGrid_bdry.m);
%   Step 2: regard mu as homotopy parameter (note that rho will change 
%       according to mu);
%   Step 3: for each fixed value of rho got in Step 1, solve corresponding
%       Davidenko differential equation using predictor-corrector method;
%   Step 4: save data and plot the outcome.
%
%   by Xige Yang
%==========================================================================

clear; close all; clc;
addpath(genpath('src'));

%% First, let's load data and re-index solutions
N=160;
% freal=sprintf('Data/start_gen_%d',N);
% fid=fopen(freal,'r');
load('Data/parameters.mat');
p.mu = 0;
mu_start = floor(p.mu);
sol_id = 1;             % solution id
prefix_input = sprintf('Data/refined_new_mu_%d-%d/N_%d.mat', mu_start, sol_id, N);
load(prefix_input);
solutions_rhos = [solutions;rhos];
num=size(solutions_rhos,2);


%% Controling mu parameter
mu_end = 2.0;
dmu = .01;
mus = mu_start:dmu:mu_end;
k = 1;
RES_ratio = .1;
solutions_test = zeros([size(solutions_rhos,1),length(mus)]);
solutions_test(:,1) = solutions_rhos(:,k);

%% fsolve options
fsolveOpts = optimset('fsolve');
fsolveOpts.Algorithm = 'levenberg-marquardt';
fsolveOpts.TolX = 1e-8;
fsolveOpts.Display = 'off';
fsolveOpts.MaxIterations = 200;

%% Main loop
tic;
p_tmp = p;
error_ratio = zeros(length(mus));   % ralative error
fprintf('Tracking %d-th solution (out of %d):\n',k ,num);
% changing mu
for l = 2:length(mus)
    p_tmp.mu = mus(l);
    %if mod(count,5) == 1
        fprintf('\t Now mu = %.3f\n',p_tmp.mu);
    %end
    solutions_test(:,l) = fsolve...
        (@(x)GS_RHS_brdy(x,p_tmp), solutions_test(:,l-1), fsolveOpts);
    error_ratio(l) = norm(abs(solutions_test(:,l)-solutions_test(:,l-1)))...
        /(norm(abs(solutions_test(:,l-1))));
    
    % Check if newton's method converges
%     if ratio > RES_ratio
%         warning('Path %d -- fsolve failed to converge at iter = %d! (res ratio = %.6f)\n',k, iter, ratio);
%         break;
%     end
end
time_end = toc;
fprintf('Iteration time = %.3fs\n',time_end);

%% Plot results and save data
%save(sprintf('Data/track_rho_%d.mat',N),'solutions_test');
ind = rhos>0;
figure(1);
plot(mus, solutions_test(end,:));
xlim(sort([mus(1),mus(end)]));
%ylim([0,max(max(solutions_test))*1.1]);
title('\mu-\rho parameter diagram');
xlabel('\mu');
ylabel('\rho');
%legend;
set(gca,'fontsize',16);
%savefig(sprintf('Imgs/mu_rho_%d.fig',N));

figure(2);
plot(mus, error_ratio);
xlim(sort([mus(1),mus(end)]));
%ylim([0,max(max(solutions_test))*1.1]);
title('Truncated error vs iteration step');
xlabel('\mu');
ylabel('ralative error');
%legend;
set(gca,'fontsize',16);