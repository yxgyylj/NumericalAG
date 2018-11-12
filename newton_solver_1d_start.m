%================ 1D GS solver with various initial guesses================
%   Problem:
%      The chemical reactions
%        S + A -> 3A, A -> C.   
%       
%      is modeled by a system of reaction diffusion equations:
%        0 = (D_A * Laplacian(A) + S.*A.^2 - (mu+rho)*A);
%        0 = (D_S * Laplacian(S) - S.*A.^2 + rho*(1-S));
%   
%   Parameters:
%        DA: diffusion rate of A;           rho: feeding rate
%        DS: diffusion rate of S;           mu:  killing rate
%
%   Dependents:
%       GS_RHS               -- right hand side function
%
%   by Xige Yang
%=========================================================================

close all;clear; clc;
addpath('src');

%% Set up system parameters
% model parameters
p.DA = .04;
p.DS = .004;
p.xl = 0.;
p.xr = 1.;

%% -- change these parameters for different solutions -- %%%
inputMsg = ['Which parameter set do you want to choose?\n'...
              '1 -- mitosis;  '...
              '2 -- corol-like;  '...
              '3 -- periodic;  '...
              '4 -- for bertiti test;\n'...
              'Your selection: \n'];
getPara = input(inputMsg);
switch getPara
    case 1  % for "mitosis" (which is expected to be homogeneous in 1d)
        p.rho = 0.0367;
        p.mu = 0.0649;

    case 2 % for corol-like structure (which is expected to be non-homogeneous in 1d)
        p.rho = 0.055;
        p.mu = 0.062;

    case 3 % for "periodic" (??? in 1d)
        p.rho = 0.018;
        p.mu = 0.051;
        
    case 4 % for bertiti test
        p.rho = 0.01;
        p.mu = 1.;

    otherwise
        error('Parameter set invalid!')
end

% %% Choose iteration algorithm
% inputMsg = ['Which fsolve algorithm do you want?\n'...
%               '1 -- levenberg-marquardt;  '...
%               '2 -- trust-region-dogleg;  '...
%               '3 -- trust-region;\n'...
%               'Your selection: \n'];
% getAlgor = input(inputMsg);
% switch getAlgor
%     case 1
%         Algor = 'levenberg-marquardt';
% 
%     case 2
%         Algor = 'trust-region-dogleg';
% 
%     case 3
%         Algor = 'trust-region';
% 
%     otherwise
%         error('Invalid algorithm for fsolve!')
% end
          
%% grid initialization and initial guess
N = 10;
xx = linspace(p.xl, p.xr, N+1)';
h = xx(2) - xx(1);
p.DA = p.DA * h^2;  % rescale DA and DS
p.DS = p.DS * h^2;

% solver options
p1.fsolveOpts = optimset('fsolve');
p1.fsolveOpts.Algorithm = 'levenberg-marquardt';
p1.fsolveOpts.TolX = 1e-8;
p1.fsolveOpts.MaxIterations = 1000;
p1.fsolveOpts.FiniteDifferenceType = 'central';
%p.fsolveOpts.Display = 'iter';

p2.fsolveOpts = optimset('fsolve');
p2.fsolveOpts.Algorithm = 'trust-region-dogleg';
p2.fsolveOpts.TolX = 1e-8;
p2.fsolveOpts.MaxIterations = 1000;
p2.fsolveOpts.FiniteDifferenceType = 'central';

p3.fsolveOpts = optimset('fsolve');
p3.fsolveOpts.Algorithm = 'trust-region';
p3.fsolveOpts.TolX = 1e-8;
p3.fsolveOpts.MaxIterations = 1000;
p3.fsolveOpts.FiniteDifferenceType = 'central';

%% Run from multiple initial guesses.
nIni_heights = 16;
nIni = nIni_heights*N;   % # of initial guesses
tic;
%solution = zeros(2*N+2,nIni);
solution_lev = [];
solution_trust = [];
solution_trust_dogleg = [];
guess = zeros(2*N+2,nIni);
TOL = 1e-4;
for k = 1:nIni
    guessA = normpdf(xx, k*h,.04*2^(mod(nIni,nIni_heights)));
    guessS = normpdf(xx, k*h,.04*2^(mod(nIni,nIni_heights)));
    guess(:,k) = [guessA; guessS];
    
    % filter solutions
    % for levenberg-marquardt method
    tmp = fsolve(@(x)GS_RHS(x,p),guess(:,k),p1.fsolveOpts);
    if k == 1
        solution_lev = tmp;
    elseif min(max(abs(tmp - solution_lev))) > TOL
        solution_lev = [solution_lev, tmp];
    end
    
    % for trust-region-dogleg method
    tmp = fsolve(@(x)GS_RHS(x,p),guess(:,k),p2.fsolveOpts);
    if k == 1
        solution_trust = tmp;
    elseif min(max(abs(tmp - solution_trust))) > TOL
        solution_trust = [solution_trust, tmp];
    end
    
    % for trust-region method
    tmp = fsolve(@(x)GS_RHS(x,p),guess(:,k),p3.fsolveOpts);
    if k == 1
        solution_trust_dogleg = tmp;
    elseif min(max(abs(tmp - solution_trust_dogleg))) > TOL
        solution_trust_dogleg = [solution_trust_dogleg, tmp];
    end
    %solution(:,k) = fsolve(@(x)GS_RHS(x,p),guess(:,k),p.fsolveOpts);
    %solution = reshape(solution,2,[]);
end
iter_end_time = toc;
fprintf('Iteration ends in %f seconds\n', iter_end_time);

%% Plot
h1 = figure;
subplot(2,2,1)
plot(xx, guess(1:N+1,:));
title('Initial guesses of A');
subplot(2,2,2)
plot(xx, guess(N+2:end,:));
title('Initial guesses of S');
subplot(2,2,3)
plot(xx, solution_lev(1:N+1,:));
title('Solutions of A');
subplot(2,2,4)
plot(xx, solution_lev(N+2:end,:));
title('Solutions of S');
%saveas(fig,'Solution.png');

h2 = figure;
subplot(3,2,1)
plot(xx, solution_lev(1:N+1,:));
title('Levenberg-Marquardt method');
subplot(3,2,2)
plot(xx, solution_lev(N+2:end,:));
title(sprintf('Number of independent steady states = %d\n', size(solution_lev,2)));
subplot(3,2,3)
plot(xx, solution_trust(1:N+1,:));
title('Trust-region method');
subplot(3,2,4)
plot(xx, solution_trust(N+2:end,:));
title(sprintf('Number of independent steady states = %d\n', size(solution_trust,2)));
subplot(3,2,5)
plot(xx, solution_trust_dogleg(1:N+1,:));
title('Trust-region-Dogleg  method');
subplot(3,2,6)
plot(xx, solution_trust_dogleg(N+2:end,:));
title(sprintf('Number of independent steady states = %d\n', size(solution_trust_dogleg,2)));
fprintf('Number of independent steady states = %d\n', size(solution_lev,2));