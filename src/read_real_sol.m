function [x, num] =read_real_sol(freal,N)
fid=fopen(freal,'r');
num=fscanf(fid,'%e',1);
x=[];
for i=1:num
    xtmp=[];
    for j=1:2*(N+1)
        xtmp(j)=fscanf(fid,'%e',1);
        tmp=fscanf(fid,'%e',1);
    end
    xtmp=xtmp';
    is_in_x=0;
    for j=1:size(x,2)
        if norm(x(:,j)-xtmp)<1e-10
            is_in_x=1;
            break
        end
    end
    if ~is_in_x
        x=[x xtmp];
    end
end
fclose(fid);
