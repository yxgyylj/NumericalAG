function [x, rho, num] = read_real_sol_bdry(freal,N,var)
if nargin == 2
    var = 0;
end

fid = fopen(freal,'r');
num = fscanf(fid,'%e',1);
x = [];
rho = [];
for i = 1:num
    xtmp = zeros(4*(N+1),1);
    for j = 1:4*(N+1)
        xtmp(j) = fscanf(fid,'%e',1);
        tmp = fscanf(fid,'%e',1);
    end
    rhotmp = fscanf(fid,'%e',1);
    if var == 1
        tmp = fscanf(fid,'%e',3);
    else
        tmp = fscanf(fid,'%e',1);
    end
    
    is_in_x = 0;
    for j = 1:size(x,2)
        if norm(x(:,j)-xtmp) < 1e-10
            is_in_x = 1;
            break
        end
    end
    if ~is_in_x
        x = [x xtmp];
        rho = [rho rhotmp];
    end
end
fclose(fid);
