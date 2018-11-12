%================= Tracking rho with parallel algorithm ===================
%
%   Step 1: load data (which is got by fineGrid_bdry.m);
%   Step 2: regard mu as homotopy parameter (note that rho will change 
%       according to mu);
%   Step 3: for each fixed value of rho got in Step 1, solve corresponding
%       Davidenko differential equation using predictor-corrector method.
%       Note that since each initial value of rho is independent from
%       others, that's where parfor loop come into play;
%   Step 4: save data and plot the outcome

%   by Xige Yang
%=========================================================================
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
mu_end = 2;
dmu = .01;
mus = mu_start:dmu:mu_end;
% solutions_track = zeros([length(mus),size(solutions_rhos)]);
% solutions_track(1,:,:) = solutions_rhos;
rhos_track = zeros(size(solutions_rhos,2),length(mus));

%% fsolve options
fsolveOpts = optimset('fsolve');
%fsolveOpts.Algorithm = 'levenberg-marquardt';
fsolveOpts.TolX = 1e-8;
fsolveOpts.Display = 'off';
fsolveOpts.MaxIterations = 200;

%% Main loop
tic;
if isempty(gcp('nocreate'))
    parpool;
end
parfor (k = 1:num,num)
    % create temperal variables for parfor
    solutions_tmp = zeros(size(solutions_rhos,1),length(mus));
    solutions_tmp(:,1) = solutions_rhos(:,k);
    p_tmp = p;
    count = 1;
    fprintf('Tracking %d-th solution (out of %d):\n',k ,num);
    % changing mu
    for mu = mus(2:end)
        if mod(count,5) == 1
            fprintf('\t Now mu = %.3f (from track %d)\n',mu, k);
        end
        p_tmp.mu = mu;
        solutions_tmp(:,count+1) = fsolve...
            (@(x)GS_RHS_brdy(x,p_tmp), solutions_tmp(:,count), fsolveOpts);
        count = count + 1;
    end
    rhos_track(k,:) = solutions_tmp(end,:);
end
time_end = toc;
fprintf('Iteration time = %.3fs\n',time_end);

%% Plot results and save data
% indicate save file paths
prefix_data = sprintf('Data/parallel/refined_new_mu_%d-%d', mu_start, sol_id);
if ~exist(prefix_data)
    mkdir(prefix_data);
end

prefix_img = sprintf('Imgs/parallel/refined_new_mu_%d-%d', mu_start, sol_id);
if ~exist(prefix_img)
    mkdir(prefix_img);
end

save(sprintf('%s/track_rho_%d_%d.mat',prefix_data, mu_start, mu_end),'rhos_track');
%ind = rhos>0 & rhos<4;
plot(mus, rhos_track);
xlim(sort([mus(1),mus(end)]));
%ylim([0,max(max(rhos_track(ind,:)))*1.1]);
% xlim([.9,1]);
% ylim([-.2,2]);
title('\mu-\rho parameter diagram');
xlabel('\mu');
ylabel('\rho');

% legend properies
str = [];
for k = 1:num
    str = [str; sprintf('Path %d',k)];
end
legend(str);
set(gca,'fontsize',24);
savefig(sprintf('%s/track_rho_%d_%d.fig',prefix_img, mu_start, mu_end));