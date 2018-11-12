%==================== RHS function for 1d Gray-Scott =====================
%   essentially all the same as GS_RHS, but with time evolution

function out = GS_RHS_t( t,x,p )

    %-- Note: Ny is usful only for 2d cases --
    x = x(:);     % delete this if in 2d case
    [Nx2, Ny] = size(x);
    h = (p.xr - p.xl)/Nx2;
    if mod(Nx2, 2) ~= 0
        warning('First dim of input matrix should be multiple of 2');
    end
    
    Nx = floor(Nx2/2);
    A = x(1:Nx);
    S = x(Nx+1:end);
    out = zeros(Nx2,1);
    
    %% right hand side functions
    
    % keep in mind -- 1:Nx represent A equation, Nx+1:Nx2 represent S equation
    % interior points
    for j = 2:Nx-1
        out(j) = (p.DA*(A(j-1) - 2*A(j) + A(j+1)))/h^2;        % diffusion
        out(j) = out(j) + A(j)^2*S(j) - (p.mu + p.rho)*A(j);   % reaction
        out(j+Nx) = (p.DS*(S(j-1) - 2*S(j) + S(j+1)))/h^2;      % diffusion
        out(j+Nx) = out(j+Nx) - A(j)^2*S(j) + p.rho*(1 - S(j));  % reaction
    end
    
    % boundary conditions
    out(1) = 3*A(1) - 4*A(2) + A(3);
    out(Nx) = 3*A(Nx) - 4*A(Nx-1) + A(Nx-2);
    out(Nx+1) = 3*S(1) - 4*S(2) + S(3);
    out(end) = 3*S(Nx) - 4*S(Nx-1) + S(Nx-2);
    
end