function GS_input_hom(N,filename,userhom,p)
    %h=0.0879645928/N;

    DA = p.DA;  DS = p.DS;    
    mu = p.mu;  rho = p.rho;
    a = p.xl;    b = p.xr;
    h=(b-a)/N;

    if exist(filename)
        delete(filename);
    end
    fid=fopen(filename,'wt');
    fprintf(fid,'CONFIG\n');

    fprintf(fid,'MPTYPE:0;\n');
    fprintf(fid,'ODEPREDICTOR: 0;\n');
    fprintf(fid,'OUTPUTLEVEL:0;\n');
    if userhom
        fprintf(fid,'USERHOMOTOPY:1;\n');
    end

    fprintf(fid,'TRACKTOLBEFOREEG: 1e-6;\n');
    fprintf(fid,'TRACKTOLDURINGEG: 1e-6;\n');
    fprintf(fid,'FINALTOL: 1e-12;\n');

    fprintf(fid,'SLICETOLBEFOREEG: 1e-6;\n');
    fprintf(fid,'SLICETOLDURINGEG: 1e-6;\n');
    fprintf(fid,'SLICEFINALTOL: 1e-12;\n');

    fprintf(fid,'SECURITYMAXNORM: 1e5;\n');
    fprintf(fid,'PathTruncationThreshold:1e5;\n');
    fprintf(fid,'EndpointFiniteThreshold:1e5;\n');
    fprintf(fid,'MINSTEPSIZEBEFOREEG:1e-5;\n');
    fprintf(fid,'MINSTEPSIZEDURINGEG:1e-5;\n');
    fprintf(fid,'NBHDRADIUS:1e-5;\n');
    fprintf(fid,'CONDNUMTHRESHOLD: 1e12;\n');
    fprintf(fid,'END;\n');
    fprintf(fid,'INPUT\n');

    ff=['function '];

    vA=[];vS=[];
    for i=1:N+1
        if i==N+1
            ff=[ff 'fA',num2str(i),','];
            ff=[ff 'fS',num2str(i),','];

        else
            ff=[ff 'fA',num2str(i),','];
            ff=[ff 'fS',num2str(i),','];
            %         if i>1
            %         ff=[ff 'fSA',num2str(i),','];
            %         end
        end
        vA=[vA;sym(['A',num2str(i)])];
        vS=[vS;sym(['S',num2str(i)])];
    end
    vxir=[];vxii=[];
    for i=1:2*(N-1)
        vxir=[vxir;sym(['xir',num2str(i)])];
        vxii=[vxii;sym(['xii',num2str(i)])];
        ff=[ff 'dfr',num2str(i),','];
        ff=[ff 'dfi',num2str(i),','];
    end
    ff=[ff 'linearr,lineari'];
    %cc=['variable_group  D,rho,mu'];
    if userhom
        variable=['variable '];
    else
        variable=['variable_group '];
    end
    for i=1:N+1
        variable=[variable,char(vA(i)),','];
    end
    for i=1:N
        variable=[variable,char(vS(i)),','];
    end
    if userhom
        variable=[variable,char(vS(N+1)),','];
    else
        variable=[variable,char(vS(N+1))];
    end
    if ~userhom
        fprintf(fid,[variable,',rho,lambda;\n']);
    else
        variable=[variable,'rho,lambda'];
    end

    if ~userhom
        variable=['variable_group '];
    end
    for i=1:2*(N-1)
        variable=[variable,char(vxir(i)),','];
    end
    for i=1:2*(N-1)-1
        variable=[variable,char(vxii(i)),','];
    end
    variable=[variable,char(vxii(2*(N-1)))];

    fprintf(fid,[variable,';\n']);


    fprintf(fid,[ff,';\n']);
    fprintf(fid,['constant h;\n']);

    if userhom
        fprintf(fid,'parameter ss;\n');
        fprintf(fid,'pathvariable tt;\n');
        fprintf(fid,'ss = tt;\n');
    end

    ftmp=['h = ',num2str(h,16),';\n'];
    fprintf(fid,ftmp);

    %fprintf(fid,['rho = ',num2str(rho,16),';\n']);
    if userhom==1
        fprintf(fid,['mu = ',num2str(rho,16),'*tt+(1-tt)*',num2str(mu,16),';\n']);
    end
    % fprintf(fid,['D = 0.1*ss+(1-ss)*',num2str(D,16), ';\n\n']);
    if userhom==2 || userhom==0
        fprintf(fid,['mu = ',num2str(mu,16),';\n']);
    end
    fprintf(fid,['DA = ',num2str(DA,16), ',\n']);
    fprintf(fid,['DS = ',num2str(DS,16), ';\n\n']);


    for i=1:N+1

        if i==1
            tmp=-3*vA(i)+4*vA(i+1)-vA(i+2);
            ftmp=['fA',num2str(i), ' = ',char(tmp),';\n'];
            fprintf(fid,ftmp);
            tmp=-3*vS(i)+4*vS(i+1)-vS(i+2);
            ftmp=['fS',num2str(i), ' = ',char(tmp),';\n'];
            fprintf(fid,ftmp);

        elseif i==N+1
            tmp=-3*vA(i)+4*vA(i-1)-vA(i-2);
            ftmp=['fA',num2str(i), ' = ',char(tmp),';\n'];
            fprintf(fid,ftmp);
            tmp=-3*vS(i)+4*vS(i-1)-vS(i-2);
            ftmp=['fS',num2str(i), ' = ',char(tmp),';\n'];
            fprintf(fid,ftmp);

        else
            tmp=sym('DA')*(vA(i+1)+vA(i-1)-2*vA(i))/sym('h')^2-vA(i)+sym('rho');
            tmp=(vS(i+1)+vS(i-1)-2*vS(i))/sym('h')^2+sym('mu')*(1+tmp);

            ftmp=['fA',num2str(i), ' = ',char(tmp),';\n'];
            fprintf(fid,ftmp);

            tmp=(vS(i+1)+vS(i-1)-2*vS(i))/sym('h')^2+sym('mu')*(1-vS(i)*vA(i)^2);
            ftmp=['fS',num2str(i), ' = ',char(tmp),';\n'];

            fprintf(fid,ftmp);

        end

    end


    A=-2*eye(N-1) + diag(ones(N-2,1),1) + diag(ones(N-2,1),-1);
    A(1)=-2/3;
    A(1,2)=2/3;
    A(end,end)=-2/3;
    A(end,end-1)=2/3;
    A=A/sym('h')^2;
    B=diag(2*vA(2:N).*vS(2:N));
    C=diag(vA(2:N).^2);

    AA=[sym('mu')*(sym('DA')*A-eye(N-1)) A;-sym('mu')*B A-sym('mu')*C];
    BB=[sym('DS')*eye(N-1)*sym('mu') eye(N-1);zeros(N-1) eye(N-1)];
    tmpr=AA*vxir+sym('lambda')*BB*vxii;
    tmpi=AA*vxii-sym('lambda')*BB*vxir;

    for i=1:2*(N-1)
        ftmp=['dfr',num2str(i), ' = ',char(tmpr(i)),';\n'];
        fprintf(fid,ftmp);
        ftmp=['dfi',num2str(i), ' = ',char(tmpi(i)),';\n'];
        fprintf(fid,ftmp);

    end

    ftmp=['linearr=',char(vxir(1)),'-1;\n'];
        fprintf(fid,ftmp);

    ftmp=['lineari=',char(vxii(1)),';\n'];
        fprintf(fid,ftmp);

    fprintf(fid,'END;');
    fclose(fid);
end
