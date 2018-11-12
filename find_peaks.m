clear; close all; clc;
addpath(genpath('src'));

%% Load data
load('Data/parameters_sep.mat','p','x');
p.N = 5;
p.H = (p.xr - p.xl)/p.N;
Ngrid = 320;
prefix_data = sprintf('Data/pure_refine_sep_N_%d/N_%d.mat', p.N, Ngrid);
load(prefix_data);

nSol = size(solutions,2);
A = solutions(1:Ngrid+1,:);
TOL = 1e-2;
ind = [];

for k = 1:nSol
    A_prime = diff(A(:,k));
    tmp = find(abs(A_prime) < TOL);
    tmp_peaks = [];
    
    % find local maxiums
    if ~isempty(tmp)
        for j = 1:length(tmp)
%             A_pp = diff(A_prime);
%             A_pp = [A_pp(1); A_pp; A_pp(end)];
%            if A_pp(tmp(j)+1) < 0
            if tmp(j) ~= 1 && tmp(j) ~= Ngrid
                if A(tmp(j),k) > A(tmp(j)+1,k) ...
                        && A(tmp(j),k) > A(tmp(j)-1,k)
                    tmp_peaks = [tmp_peaks, tmp(j)];
                end
            end
        end
    end
    Npeaks = length(tmp_peaks);
    
    % find solutions with 4 peaks
    if Npeaks == 4
        ind = [ind, k];
    end
end

figure('Position',[300,300,400,400])
if ~isempty(ind)
    for j = 0:3
        subplot(2,2,j+1)
        plot(finer_grid,A(:,ind(4*j+1:4*j+4)),'linewidth',2);
        xlabel('x');
        title('A');
    end
else
    warning('No solutions with 4 peaks!');
end
        