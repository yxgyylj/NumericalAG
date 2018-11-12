%========================= Make a finer grid =============================
%   Coarse grid size is H = (xr - xl)/N
%   Finer grid sizes are h = H/2^j for some non-neg integer j
%
%       Step 1: get the solutions on coarse grid;
%       Step 2: interpolate the solution onto a finer grid;
%       Step 3: use the interpolated solution as initial guess, use 
%               Newton's method to get new solutions;
%       Step 4: go back to Step 2 for a preallocated step numbers.
%
%   by Xige Yang
%=========================================================================

clear; close all; clc;
addpath(genpath('src'));

%% Parameters
% p.N = 5;                % # of subdomains
% p.DA =0.1/p.N^2;        % diffusion of A
% p.DS = 1./p.N^2;        % diffusion of S
% p.mu=1.;                % killing rate
% p.rho = 0.01;           % feeding rate
% p.xl = 0.;              % left boundary
% p.xr = 1.;              % right boundary
% p.H = (p.xr - p.xl)/p.N;
% x=linspace(p.xl,p.xr,p.N+1)';
load('Data/parameters_sep.mat','p','x');
p.N = 5;
p.H = (p.xr - p.xl)/p.N;


%% Extract solutions
freal=sprintf('bertiniFile/Sep/input_%d/real_finite_solutions',p.N);
[init_sol, nSol] = read_real_sol(freal,p.N);
% fprintf('Number of solutions is: %d\n', size(SA(1,:),2));

%% Save file paths
prefix_data = sprintf('Data/pure_refine_sep_N_%d', p.N);
if ~exist(prefix_data)
    mkdir(prefix_data);
end

prefix_img = sprintf('Imgs/pure_refine_sep_N_%d', p.N);
if ~exist(prefix_img)
    mkdir(prefix_img);
end

%% Get finer grids
% initialize
h = p.H;
solutions = init_sol;
nSol_finer = nSol;

% start main loop
tic;
for j = 1:6  
    % Plot current solution
    xx = (p.xl:h:p.xr)';
    nGrid_old = length(xx);
    figure(1);
    subplot(2,1,1)
    plot(xx,solutions(1:nGrid_old,:))
    xlabel('x');
    ylabel('A');
    title(sprintf('Solutions in N=%d subdomains', p.N*2^(j-1)));
    set(gca,'fontsize',20);
    subplot(2,1,2)
    plot(xx,solutions(nGrid_old+1:end,:))
    xlabel('x');
    ylabel('S');
    title(sprintf('Number of solutions = %d', nSol_finer));
    set(gca,'fontsize',20);
    savefig(sprintf('%s/N_%d.fig', prefix_img, p.N*2^(j-1)));

    % drop negative solutions
    k = 1;
    large_sol = 2e3;
    while k <= size(solutions,2)
        if min(solutions(1:2*nGrid_old,k)) < -.1
            solutions(:,k) = [];
        else
            k = k + 1;
        end
    end    
    
    % iteration vectors and variables onto a finer grid
    h = h/2;
    finer_grid = p.xl:h:p.xr;
    nGrid = length(finer_grid);
    
    %% Once we get the finer grid, iterate along all previous solutions
    [row_sol, col_sol] = size(solutions);
    solutions = [solutions; zeros(row_sol-2, col_sol)];
    sol_tmp = [];
    
    % update initial guesses
    for k = 1:col_sol
        tmp_sol_A = solutions(1:nGrid_old,k);
        tmp_sol_S = solutions(nGrid_old+1:2*nGrid_old,k);
        guess_A = spline(xx, tmp_sol_A, finer_grid); 
        guess_S = spline(xx, tmp_sol_S, finer_grid); 
        guess = [guess_A(:); guess_S(:)];
        sol_refine = myNewton(@GS_RHS, guess, p);
        TOL = 1e-6;
        
        % filter solutions
        if k == 1
            sol_tmp = sol_refine;
        elseif min(max(abs(sol_tmp - sol_refine))) > TOL
            sol_tmp = [sol_tmp, sol_refine];
        end
    end
    
        %solution(:,k) = refine_sol;
    
    % update refined solutions
    solutions = sol_tmp;
    nGrid_old = nGrid;
    nSol_finer = size(solutions,2);
    save(sprintf('Data/solutions_%d.mat',p.N*2^j),'solutions','finer_grid');
    
    % generate start files
    MyFilename = sprintf('Data/start_gen_%d',p.N*2^j);
    bertini_start(MyFilename,solutions);
end

end_time = toc;
fprintf('Elapsed time = %f seconds \n', end_time);
save(sprintf('%s/N_%d.mat', prefix_data, p.N*2^j),'solutions','finer_grid');

%% Plot final outcome
figure(1);
subplot(2,1,1)
plot(finer_grid,solutions(1:nGrid,:))
xlabel('x');
ylabel('A');
title(sprintf('Solutions in N=%d subdomains', p.N*2^j));
set(gca,'fontsize',16);
subplot(2,1,2)
plot(finer_grid,solutions(nGrid+1:2*nGrid,:))
xlabel('x');
ylabel('S');
title(sprintf('Number of solutions = %d', nSol_finer));
set(gca,'fontsize',16);
savefig(sprintf('%s/N_%d.fig', prefix_img, p.N*2^j));