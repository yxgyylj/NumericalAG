function GS_input_lin(N,filename,parahom,p)
    % problem parameters
    DA = p.DA;  DS = p.DS;    
    mu = p.mu;  rho = p.rho;
    a = p.xl;   b = p.xr;
    h=(b-a)/N;
    rho_start = .01;
    
    if exist(filename)
        delete(filename);
    end
    fid=fopen(filename,'wt');
    
    % configuring
    fprintf(fid,'CONFIG\n');

    fprintf(fid,'MPTYPE:0;\n');
    fprintf(fid,'ODEPREDICTOR: 0;\n');
    fprintf(fid,'OUTPUTLEVEL:0;\n');
    if parahom
        fprintf(fid,'PARAMETERHOMOTOPY:2;\n');
    end

    fprintf(fid,'TRACKTOLBEFOREEG: 1e-6;\n');
    fprintf(fid,'TRACKTOLDURINGEG: 1e-6;\n');
    fprintf(fid,'FINALTOL: 1e-12;\n');

    fprintf(fid,'SECURITYMAXNORM: 1e5;\n');
    fprintf(fid,'PathTruncationThreshold:1e5;\n');
    fprintf(fid,'EndpointFiniteThreshold:1e5;\n');
    fprintf(fid,'MINSTEPSIZEBEFOREEG:1e-5;\n');
    fprintf(fid,'MINSTEPSIZEDURINGEG:1e-5;\n');
%    fprintf(fid,'NBHDRADIUS:1e-5;\n');
    fprintf(fid,'CONDNUMTHRESHOLD: 1e12;\n');
    fprintf(fid,'END;\n');
    fprintf(fid,'INPUT\n');

    ff=['function '];

    vA=[];vS=[];
    vAl=[];vSl=[];
    for k=1:N+1
        if k==N+1
            ff=[ff 'fA',num2str(k),','];
            ff=[ff 'fS',num2str(k),','];
            ff=[ff 'fAl',num2str(k),','];
            ff=[ff 'fSl',num2str(k),','];

        else
            ff=[ff 'fA',num2str(k),','];
            ff=[ff 'fS',num2str(k),','];
            ff=[ff 'fAl',num2str(k),','];
            ff=[ff 'fSl',num2str(k),','];
            %         if i>1
            %         ff=[ff 'fSA',num2str(i),','];
            %         end
        end
        vA=[vA;sym(['A',num2str(k)])];
        vS=[vS;sym(['S',num2str(k)])];
        vAl=[vAl;sym(['Al',num2str(k)])];
        vSl=[vSl;sym(['Sl',num2str(k)])];
    end
    ff=[ff 'linear'];
    %cc=['variable_group  DA,DB,rho,mu'];
    if parahom
        variable=['variable_group '];
    else
        variable=['variable_group '];
    end
    for k=1:N+1
        variable=[variable,char(vA(k)),','];
    end
    for k=1:N
        variable=[variable,char(vS(k)),','];
    end
    variable=[variable,char(vS(N+1)),','];

%     if parahom
%         variable=['variable_group '];
%     end
    for k=1:N+1
        variable=[variable,char(vAl(k)),','];
    end
    for k=1:N
        variable=[variable,char(vSl(k)),','];
    end
    fprintf(fid,[variable,char(vSl(N+1)),';\n']);



    % variable=['variable_group '];
    % for i=1:N-2
    %     variable=[variable,char(vSA(i)),','];
    % end
    % variable=[variable,char(vSA(N-1))];
    % fprintf(fid,[variable,';\n']);


    fprintf(fid,[ff,';\n']);
    fprintf(fid,'constant h,mu,DA,DS,');

    if parahom
        fprintf(fid,'parameter rho;\n');
%        fprintf(fid,'pathvariable tt;\n');
    else
        fprintf(fid,'rho;\n');
    end
 

    % define variable groups
    % homotopy tracking used

    ftmp=['h = ',num2str(h,16),';\n'];
    fprintf(fid,ftmp);

    %fprintf(fid,['rho = ',num2str(rho,16),';\n']);
%     if parahom==1
%         fprintf(fid,['rho = ',num2str(rho,16),'*(1-ss)+ss*',num2str(rho_start,16),';\n']);
%     end
%     % fprintf(fid,['rho = 0.1*ss+(1-ss)*',num2str(rho,16), ';\n\n']);
%     if parahom==2 || parahom==0
%         fprintf(fid,['mu = ',num2str(mu,16),';\n']);
%     end
    fprintf(fid,['mu = ',num2str(mu,16),';\n']);
    fprintf(fid,['DA = ',num2str(DA,16), ';\n']);
    fprintf(fid,['DS = ',num2str(DS,16), ';\n\n']);


    for k=1:N+1
        fprintf(fid,'\n%% The %d-th grid point \n',k);
        % left boundary
        if k==1
            tmp=-3*vA(k)+4*vA(k+1)-vA(k+2);
            ftmp=['fA',num2str(k), ' = ',char(tmp),';\n'];
            fprintf(fid,ftmp);
            tmp=-3*vS(k)+4*vS(k+1)-vS(k+2);
            ftmp=['fS',num2str(k), ' = ',char(tmp),';\n'];
            fprintf(fid,ftmp);
            
            tmp=-3*vAl(k)+4*vAl(k+1)-vAl(k+2);
            ftmp=['fAl',num2str(k), ' = ',char(tmp),';\n'];
            fprintf(fid,ftmp);
            tmp=-3*vSl(k)+4*vSl(k+1)-vSl(k+2);
            ftmp=['fSl',num2str(k), ' = ',char(tmp),';\n'];
            fprintf(fid,ftmp);

        % right boundary
        elseif k==N+1
            tmp=-3*vA(k)+4*vA(k-1)-vA(k-2);
            ftmp=['fA',num2str(k), ' = ',char(tmp),';\n'];
            fprintf(fid,ftmp);
            tmp=-3*vS(k)+4*vS(k-1)-vS(k-2);
            ftmp=['fS',num2str(k), ' = ',char(tmp),';\n'];
            fprintf(fid,ftmp);
            
            tmp=-3*vAl(k)+4*vAl(k-1)-vAl(k-2);
            ftmp=['fAl',num2str(k), ' = ',char(tmp),';\n'];
            fprintf(fid,ftmp);
            tmp=-3*vSl(k)+4*vSl(k-1)-vSl(k-2);
            ftmp=['fSl',num2str(k), ' = ',char(tmp),';\n'];
            fprintf(fid,ftmp);


        else
            % interior A and its linearized equations
            tmp = sym('DA')*(vA(k+1)+vA(k-1)-2*vA(k))/sym('h')^2;
            tmp = tmp+vS(k)*vA(k)^2-(sym('rho')+sym('mu'))*vA(k);
            ftmp = ['fA',num2str(k), ' = ',char(tmp),';\n'];
            fprintf(fid,ftmp);
            
%             tmp=sym('DA')*(vA(k+1)+vA(k-1)-2*vA(k))/sym('h')^2-vA(k)+sym('rho');
%             tmp=(vS(k+1)+vS(k-1)-2*vS(k))/sym('h')^2+sym('mu')*(1+tmp);
%             ftmp=['fA',num2str(k), ' = ',char(tmp),';\n'];
%             fprintf(fid,ftmp);
%             tmp=(vS(k+1)+vS(k-1)-2*vS(k))/sym('h')^2+sym('mu')*(1-vS(k)*vA(k)^2);
%             ftmp=['fS',num2str(k), ' = ',char(tmp),';\n'];
%             fprintf(fid,ftmp);

            % interior S and its linearized equations
            tmp = sym('DS')*(vS(k+1)+vS(k-1)-2*vS(k))/sym('h')^2;
            tmp = tmp + sym('DA')*(vA(k+1)+vA(k-1)-2*vA(k))/sym('h')^2;
            tmp = tmp - (sym('rho')+sym('mu'))*vA(k) - sym('rho')*(1-vS(k));
            ftmp = ['fS',num2str(k), ' = ',char(tmp),';\n'];
            fprintf(fid,ftmp);
            
            tmp=(sym('DA')*(vAl(k+1)+vAl(k-1)-2*vAl(k))/sym('h')^2)...
                +sym('DS')*(vSl(k+1)+vSl(k-1)-2*vSl(k))/sym('h')^2;
            tmp=tmp-(sym('mu')+sym('rho'))*vAl(k)-sym('rho')*vSl(k);
            ftmp=['fAl',num2str(k), ' = ',char(tmp),';\n'];
            fprintf(fid,ftmp);
            tmp=sym('DS')*(vSl(k+1)+vSl(k-1)-2*vSl(k))/sym('h')^2;
            tmp=tmp-vSl(k)*vA(k)^2-2*vA(k)*vS(k)*vAl(k)-sym('rho')*vSl(k);
            ftmp=['fSl',num2str(k), ' = ',char(tmp),';\n'];
            fprintf(fid,ftmp);

        end

    end
%     ftmp=['linear = ',char(vpa([vAl;vSl].'*0.5/(N+1)*ones(2*(N+1),1))),'-1;\n'];
      tmp=0;
      for i=1:N
         tmp=tmp+((vAl(i)^2+vAl(i+1)^2)/2+(vSl(i)^2+vSl(i+1)^2)/2)*sym('h');
      end
    
%    ftmp=['linear = ',char(tmp),'-1;\n'];
    ftmp=['linear=',char(vAl(2)),'-1;\n'];
    fprintf(fid,ftmp);

    fprintf(fid,'END;');
    fclose(fid);
end
