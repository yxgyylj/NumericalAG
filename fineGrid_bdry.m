%================= Make a finer grid for parameter boundaries==============
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
p.N = 5;                % # of subdomains
p.DA =0.1/p.N^2;        % diffusion of A
p.DS = 1./p.N^2;        % diffusion of S
p.mu = 0.;                % killing rate
p.rho = 0.01;           % feeding rate
p.xl = 0.;              % left boundary
p.xr = 1.;              % right boundary
p.H = (p.xr - p.xl)/p.N;
x=linspace(p.xl,p.xr,p.N+1)';
mu_start = floor(p.mu);
sol_id = 3;             % solution id
save('Data/parameters.mat','p','x', 'sol_id');


%% Extract solutions
prefix = 'bertiniFile/input_bdry_5/new_solutions_mu';
freal = sprintf('%s_%d-%d',prefix,mu_start,sol_id);
[init_sol, rho_sin, nSol] = read_real_sol_bdry(freal,p.N);

%% Save file paths
prefix_data = sprintf('Data/refined_new_mu_%d-%d', mu_start, sol_id);
if ~exist(prefix_data)
    mkdir(prefix_data);
end

prefix_img = sprintf('Imgs/refined_new_mu_%d-%d', mu_start, sol_id);
if ~exist(prefix_img)
    mkdir(prefix_img);
end
% fprintf('Number of solutions is: %d\n', size(SA(1,:),2));

%% Get finer grids
% initialize
h = p.H;
solutions = init_sol;
rhos = rho_sin;
nSol_finer = nSol;

% start main loop
tic;
for j = 1:6  
    % Plot current solution
    xx = (p.xl:h:p.xr)';
    nGrid_old = length(xx);
    
    figure(1);
    set(gca,'fontsize',16);
    subplot(2,2,1)
    plot(xx,solutions(1:nGrid_old,:))
    xlabel('x');
    ylabel('A');
    title(sprintf('Solutions in N=%d subdomains', p.N*2^(j-1)));
    set(gca,'fontsize',16);
    
    subplot(2,2,3)
    plot(xx,solutions(nGrid_old+1:2*nGrid_old,:))
    xlabel('x');
    ylabel('S');
    title(sprintf('Number of solutions = %d', nSol_finer));
    set(gca,'fontsize',16);
    
    subplot(2,2,2)
    plot(xx,solutions(2*nGrid_old+1:3*nGrid_old,:))
    xlabel('x');
    ylabel('A1');
    title('First order terms');
    set(gca,'fontsize',16);
    
    subplot(2,2,4)
    plot(xx,solutions(3*nGrid_old+1:4*nGrid_old,:))
    xlabel('x');
    ylabel('S1');
    set(gca,'fontsize',16);
    fprintf('When N=%d, rhos=[%s] \n', p.N*2^(j-1),num2str(rhos));
    savefig(sprintf('%s/N_%d.fig', prefix_img, p.N*2^(j-1)));

    % drop large and negative solutions
    k = 1;
    large_sol = 2e3;
    large_sol_diff = 1e4;
    while k <= size(solutions,2)
        if max(abs(solutions(1:2*nGrid_old,k))) > large_sol ...
                || min(solutions(1:2*nGrid_old,k)) < -.1 ...
                || max(abs(solutions(2*nGrid_old:end-1,k))) > large_sol_diff ...
                || solutions(end,k) < -20
            solutions(:,k) = [];
        else
            k = k + 1;
        end
    end    
    nSol_finer = size(solutions,2);
    
    % iteration vectors and variables onto a finer grid
    h = h/2;
    finer_grid = p.xl:h:p.xr;
    nGrid = length(finer_grid);
    
    %% Once we get the finer grid, iterate along all previous solutions
    [row_sol, col_sol] = size(solutions);
    %solutions = [solutions; zeros(row_sol-2, col_sol)];
    sol_tmp = [];
    %rho_tmp = rhos(1);
    res = zeros(1,nSol_finer);
    
    % update initial guesses
    for k = 1:nSol_finer
        tmp_sol_A = solutions(1:nGrid_old,k);
        tmp_sol_S = solutions(nGrid_old+1:2*nGrid_old,k);
        tmp_sol_Al = solutions(2*nGrid_old+1:3*nGrid_old,k);
        tmp_sol_Sl = solutions(3*nGrid_old+1:4*nGrid_old,k);
        guess_A = spline(xx, tmp_sol_A, finer_grid); 
        guess_S = spline(xx, tmp_sol_S, finer_grid); 
        guess_Al = spline(xx, tmp_sol_Al, finer_grid); 
        guess_Sl = spline(xx, tmp_sol_Sl, finer_grid); 
        guess = [guess_A(:); guess_S(:); guess_Al(:); guess_Sl(:); rhos(k)];
        sol_refine = myNewton(@GS_RHS_brdy, guess, p);
        res(k) = norm(GS_RHS_brdy(sol_refine, p));
        TOL1 = 1e-6;
        TOL2 = 1e-2;
        
        % filtering solutions
        if k == 1
            sol_tmp = sol_refine;
        elseif min(max(abs(sol_tmp - sol_refine))) > TOL1 && res(k) < TOL2
            sol_tmp = [sol_tmp, sol_refine];
            %rho_tmp = [rho_tmp, rhos(k)];
        end
    end
    
        %solutions(:,k) = refine_sol;
    
    % update refined solutions
    solutions = sol_tmp(1:end-1,:);
    
    % delete the first solution if not converge by Newton's method
    if res(1) >= TOL2
        solutions(:,1) = [];
        res(1) = [];
    end
    
    % save data
    rhos = sol_tmp(end,:);
    nGrid_old = nGrid;
    nSol_finer = size(solutions,2);
    save(sprintf('%s/N_%d.mat', prefix_data, p.N*2^j),'solutions','rhos','finer_grid');
    %prefix = sprintf('Data/start_gen_brdy_%d',p.N*2^j);
    %bertini_start(MyFilename,solutions);
end

end_time = toc;
fprintf('Elapsed time = %f seconds \n', end_time);
sort(solutions(1:rhos));
save(sprintf('%s/N_%d.mat', prefix_data, p.N*2^j),'solutions','rhos','finer_grid');

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
savefig(sprintf('%s/N_%d.fig', prefix_img, p.N*2^j));

figure(2);
subplot(2,1,1)
plot(1:nSol_finer,solutions(end,:)','b*');
title('\rho values');
xlabel('solution id');
ylabel('rho');
set(gca,'fontsize',16);

subplot(2,1,2)
plot(1:nSol_finer,res,'r*');
title(['redial plot of Newton\prime', 's method']);
xlabel('solution id');
ylabel('L^2 norm');
set(gca,'fontsize',16);
savefig(sprintf('%s/error_plot.fig', prefix_img));