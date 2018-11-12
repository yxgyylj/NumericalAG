%============= 1D GrayScott Newton solver with various guesses ============
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
p.DA = .1;
p.DS = 1.;
p.xl = 0.;
p.xr = 1.;

%%% -- change these parameters for different solutions -- %%%
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

%% grid initialization and initial guess
N = 128;
xx = linspace(p.xl, p.xr, N+1);
h = xx(2) - xx(1);
p.DA = p.DA * h^2;  % rescale DA and DS
p.DS = p.DS * h^2;
guessA = normpdf(xx, .8,.25);
guessS = normpdf(xx, .2,.25);
guess = [guessA; guessS];

% solver options
p.fsolveOpts = optimset('fsolve');
p.fsolveOpts.Algorithm = 'levenberg-marquardt';
p.fsolveOpts.TolX = 1e-8;
%p.fsolveOpts.Display = 'iter';

%% Run from the initial guess.
tic;
solution = fsolve(@(x)GS_RHS(x,p),guess,p.fsolveOpts);
%solution = reshape(solution,2,[]);
iter_end_time = toc;

%% Plot
fig1 = figure;
subplot(2,1,1)
plot(xx, guess);
title('Initial guess');
legend('A', 'S');
subplot(2,1,2)
plot(xx, solution(1:N+1));
hold on;
plot(xx, solution(N+2:end));
title('Solution via Newton iteration');
legend('A', 'S');
%saveas(fig,'Solution.png');
fprintf('Iteration ends in %f seconds\n', iter_end_time);
