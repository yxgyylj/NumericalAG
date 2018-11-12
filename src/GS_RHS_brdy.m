%==================== RHS function for 1d Gray-Scott =====================

function out = GS_RHS_brdy(x,p)

    %-- Note: Ny is usful only for 2d cases --
    x = x(:);     % delete this if in 2d case
    [Nx4, Ny] = size(x);    % not used for 2d cases
    
    Nx = floor(Nx4/4);
    h = (p.xr - p.xl)/(Nx-1);
    A = x(1:Nx);
    S = x(Nx+1:2*Nx);
    Al = x(2*Nx+1:3*Nx);
    Sl = x(3*Nx+1:4*Nx);
    rho_var = x(end);
    
    
    if mod(Nx4, 4) == 0
        rho_var = p.rho;
    elseif mod(Nx4, 4) ~= 1
        warning('First dim of input matrix should be of form 4n or 4n+1!');
    end
    
    out = zeros(Nx4,1);
    
    %% right hand side functions
    
    % keep in mind -- 1:Nx represent A equation, Nx+1:2Nx represent S equation
    %       2Nx+1:3Nx represent Al equation, 3Nx+1:4Nx represent Sl equation
    % interior points
    for j = 2:Nx-1
        out(j) = (p.DA*(A(j-1) - 2*A(j) + A(j+1)))/h^2;        % diffusion
        out(j) = out(j) + A(j)^2*S(j) - (p.mu + rho_var)*A(j);   % reaction
        out(j+Nx) = (p.DS*(S(j-1) - 2*S(j) + S(j+1)))/h^2;      % diffusion
        out(j+Nx) = out(j+Nx) - A(j)^2*S(j) + rho_var*(1 - S(j));  % reaction
        
        % diffusion and reaction for first order perturbation terms
        out(j+2*Nx) = (p.DA*(Al(j-1) - 2*Al(j) + Al(j+1)))/h^2;      
        out(j+2*Nx) = out(j+2*Nx) + A(j)^2*Sl(j) + 2*S(j)*A(j)*Al(j) ...
             - (p.mu + rho_var)*Al(j);
        out(j+3*Nx) = (p.DS*(Sl(j-1) - 2*Sl(j) + Sl(j+1)))/h^2;
        out(j+3*Nx) = out(j+3*Nx) - A(j)^2*Sl(j) - 2*S(j)*A(j)*Al(j) - rho_var*Sl(j);
    end
    
    % boundary conditions
    out(1) = 3*A(1) - 4*A(2) + A(3);
    out(Nx) = 3*A(Nx) - 4*A(Nx-1) + A(Nx-2);
    out(Nx+1) = 3*S(1) - 4*S(2) + S(3);
    out(2*Nx) = 3*S(Nx) - 4*S(Nx-1) + S(Nx-2);
    
    out(2*Nx+1) = 3*Al(1) - 4*Al(2) + Al(3);
    out(3*Nx) = 3*Al(Nx) - 4*Al(Nx-1) + Al(Nx-2);
    out(3*Nx+1) = 3*Sl(1) - 4*Sl(2) + Sl(3);
    out(4*Nx) = 3*Sl(Nx) - 4*Sl(Nx-1) + Sl(Nx-2);
    
    if mod(Nx4, 4) == 1
        out(4*Nx+1) = Al(2) - 1;
    end
end