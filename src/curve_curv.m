function out = curve_curv(x,z)
% Curvature of an 1d curve z = z(x,y)
%   inputs:
%       x : 1d meshgrid
%       z : vector that have same length of x
    [ncol, nrow] = size(z);
    x = x(:); z = reshape(z,length(x),[]);
    zx = gradient(z)./gradient(x);
    zxx = gradient(zx)./gradient(x);
    out = zxx./(1+zx.^2);
    out = reshape(out,ncol,nrow);
end