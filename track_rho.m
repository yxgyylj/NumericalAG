%======================== Start tracking for rho =========================
%
%   Step 1: load data (which is got by fineGrid_bdry.m);
%   Step 2: regard mu as homotopy parameter (note that rho will change 
%       according to mu);
%   Step 3: for each fixed value of rho got in Step 1, solve corresponding
%       Davidenko differential equation using predictor-corrector method;
%   Step 4: save data and plot the outcome.
%
%   by Xige Yang
%=========================================================================

clear; close all; clc;
addpath(genpath('src'));

%% First, let's load data and re-index solutions
N=320;
% freal=sprintf('Data/start_gen_%d',N);
% fid=fopen(freal,'r');
load('Data/parameters.mat');
load(sprintf('Data/solution_brdy_%d.mat',N));
[rhos, ind] = sort(rhos);
solutions_rhos = [solutions(:,ind);rhos];   % solutions and rhos
num=size(solutions_rhos,2);

%% Controling mu parameter
mu_start = p.mu;
mu_end = .90;
dmu = -.01;
mus = mu_start:dmu:mu_end;
solutions_track = zeros([size(solutions_rhos),length(mus)]);
solutions_track(:,:,1) = solutions_rhos;

%% fsolve options
fsolveOpts = optimset('fsolve');
%fsolveOpts.Algorithm = 'levenberg-marquardt';
fsolveOpts.TolX = 1e-8;
fsolveOpts.Display = 'off';
fsolveOpts.MaxIterations = 100;

%% Main loop
tic;
for k = 1:num
    p_tmp = p;
    count = 1;
    fprintf('Tracking %d-th solution (out of %d):\n',k ,num);
    % changing mu
    for mu = mus(2:end)
        %if mod(count,5) == 1
            fprintf('\t Now mu = %.3f\n',mu);
        %end
        p_tmp.mu = mu;
        solutions_track(:,k,count+1) = fsolve...
            (@(x)GS_RHS_brdy(x,p_tmp), solutions_track(:,k,count), fsolveOpts);
        count = count + 1;
    end
end
time_end = toc;
fprintf('Iteration time = %.3fs\n',time_end);

%% Plot results and save data
save(sprintf('Data/track_rho_%d.mat',N),'solutions_track');
ind = rhos>0;
plot(mus, squeeze(solutions_track(end,ind,:)));
xlim(sort([mus(1),mus(end)]));
ylim([0,max(max(squeeze(solutions_track(end,ind,:))))*1.1]);
title('\mu-\rho parameter regium');
xlabel('\mu');
ylabel('\rho');
legend;
set(gca,'fontsize',16);
savefig(sprintf('Imgs/mu_rho_%d.fig',N));