% Newton's method on Gray-Scott equation

function [out, num] = myNewton(RHS, guess, p)
%% initial guess
[nRow,nCol] = size(guess);
% guessA = normpdf(xx, .8,.25);
% guessS = normpdf(xx, .2,.25);
% guess = [guessA; guessS];

% solver options
fsolveOpts = optimset('TolFun',1e-30);
%fsolveOpts.Algorithm = 'levenberg-marquardt';
fsolveOpts.TolX = 1e-8;
fsolveOpts.Display = 'off';
%fsolveOpts.MaxIterations = 200;

%% Run from the initial guess.
% only select the solutions that converge
tolerance = 1e-2;
num = 0;
for j = 1:nCol
    [out, fval, exitflag] = fsolve(@(x)RHS(x,p),guess,fsolveOpts);
%     if exitflag > 0 && j == 1
%         out = temp;
%         num = num + 1;
%         continue
%     end
    
%    res = abs(out - temp);
%     if exitflag > 0 && max(res(:)) > tolerance
%         out = [out, temp];
%         num = num + 1;
%     end
    %solution = reshape(solution,2,[]);
end

end