%=================== Jacobian of 1d Gray-Scott Problem ====================

function J = GS_Jacobian(N,x,p)
    h = (p.xl - p.xr)/N;
    A = x(1:N);
    S = x(N+1:end);
    
    JA = -2*eye(N) + diag(ones(N-1,1),1) + diag(ones(N-1,1),-1);
    JA = p.DA*JA/h^2+diag(2*A.*S-p.mu-p.rho);
    JA(1,1) = -3;   JA(1,2) = 4;    JA(1,3) = -1;
    JA(end,end) = -3;   JA(end,end-1) = 4;    JA(end,end-2) = -1;
    
    JS = -2*eye(N) + diag(ones(N-1,1),1) + diag(ones(N-1,1),-1);
    JS = p.DS*JS/h^2-p.mu*diag(A.^2);
    JS(1,1) = -3;   JS(1,2) = 4;    JS(1,3) = -1;
    JS(end,end) = -3;   JS(end,end-1) = 4;    JS(end,end-2) = -1;

    JAS=diag(A.^2);    
    JAS(1,1) = 0;   JAS(end,end) = 0;
    
    JSA=diag(-2*A.*S);    
    JSA(1,1) = 0;   JSA(end,end) = 0;

    J=[JA JAS;JSA JS];
end


