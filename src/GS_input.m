function GS_input(N,filename,Homflag,p)

DA = p.DA;  DS = p.DS;    
mu = p.mu;  rho = p.rho;
a = p.xl;    b = p.xr;

h=(b-a)/N;
fid=fopen(filename,'wt');

%% Configure part
fprintf(fid,'CONFIG\n');

fprintf(fid,'MPTYPE:2;\n');
fprintf(fid,'SecurityLevel:1;\n');
fprintf(fid,'FINALTOL: 1e-16;\n');
fprintf(fid,'AMPMaxPrec: 2048;\n');
fprintf(fid,'NBHDRADIUS:1e-5;\n');

if Homflag
    fprintf(fid,'USERHOMOTOPY:1;\n');
end

fprintf(fid,'END;\n\n');

%% Inputs
fprintf(fid,'INPUT\n');

% claim variables and functions
ff=['function '];

vA=[];
vS=[];
for j=1:N+1
    if j==N+1
        ff=[ff 'fA',num2str(j),','];
        ff=[ff 'fS',num2str(j)];
    else
        ff=[ff 'fA',num2str(j),','];
        ff=[ff 'fS',num2str(j),','];
    end
    vA=[vA;sym(['A',num2str(j)])];
    vS=[vS;sym(['S',num2str(j)])];
    if j>1 && j<N+1
%     vSA=[vSA;sym(['SA',num2str(i)])];
    end
end

cc=['constant DA, DS, rho, mu'];

% define variable groups
    if Homflag
        variable=['variable '];
    else
        variable=['variable_group '];
    end

    for j=1:N+1
        variable=[variable,char(vA(j)),','];    
    end
    
%     variable=[variable,char(vA(N+1))];    
%     fprintf(fid,[variable,';\n']);
% 
%     if Homflag
%         variable=['variable '];
%     else
%         variable=['variable_group '];
%     end
    
    for j=1:N
        variable=[variable,char(vS(j)),','];    
    end
    variable=[variable,char(vS(N+1))];
    fprintf(fid,[variable,';\n']);

fprintf(fid,[ff,';\n']);
fprintf(fid,[cc,';\n']);
if Homflag
    fprintf(fid,'parameter ss;\n');
    fprintf(fid,'pathvariable tt;\n');
    fprintf(fid,'ss = tt;\n');
end

ftmp=['h = ',num2str(h,16),';\n'];
fprintf(fid,ftmp);
fprintf(fid,['DA =',num2str(DA,16), ';\n']);
fprintf(fid,['DS =',num2str(DS,16), ';\n']);
fprintf(fid,['rho = ',num2str(rho,16),';\n']);
fprintf(fid,['mu = ',num2str(mu,16),';\n\n']);

for j=1:N+1
    
    % left boundary
    if j == 1
        tmp = -3*vA(j)+4*vA(j+1)-vA(j+2);
        ftmp = ['fA',num2str(j), ' = ',char(tmp),';\n'];
        fprintf(fid,ftmp);
        tmp = -3*vS(j)+4*vS(j+1)-vS(j+2);
        ftmp = ['fS',num2str(j), ' = ',char(tmp),';\n'];
        fprintf(fid,ftmp);
    
    % right boundary
    elseif j == N+1
        tmp = -3*vA(j)+4*vA(j-1)-vA(j-2);
        ftmp = ['fA',num2str(j), ' = ',char(tmp),';\n'];
        fprintf(fid,ftmp);
        tmp = -3*vS(j)+4*vS(j-1)-vS(j-2);
        ftmp = ['fS',num2str(j), ' = ',char(tmp),';\n'];
        fprintf(fid,ftmp);
        
    else
        % interior A equations
        %tmp=(vS(i+1)+vS(i-1)-2*vS(i))/(sym('h')^2*sym('mu'))+1;
        tmp = sym('DA')*(vA(j+1)+vA(j-1)-2*vA(j))/sym('h')^2;
        tmp = tmp+vS(j)*vA(j)^2-(sym('rho')+sym('mu'))*vA(j);
        ftmp = ['fA',num2str(j), ' = ',char(tmp),';\n'];
        fprintf(fid,ftmp);
        
        % interior S equations
        tmp = sym('DS')*(vS(j+1)+vS(j-1)-2*vS(j))/sym('h')^2;
        tmp = tmp + sym('DA')*(vA(j+1)+vA(j-1)-2*vA(j))/sym('h')^2;
        tmp = tmp - (sym('rho')+sym('mu'))*vA(j) + sym('rho')*(1-vS(j));
        ftmp = ['fS',num2str(j), ' = ',char(tmp),';\n'];
        fprintf(fid,ftmp);
        
%         tmp=vSA(i-1)-vS(i)*vA(i)^2;
%         ftmp=['fSA',num2str(i), ' = ',char(tmp),';\n'];
%         
%         fprintf(fid,ftmp);

    end
    
end
fprintf(fid,'END;');
fclose(fid);