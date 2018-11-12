%====================== Generate bertini start system ====================

function bertini_start(filename,solution)

    if exist(filename)
        delete filename;
    end
    
    lSol = size(solution,1);    % length of solution
    nSol = size(solution,2);    % solution number
    fid=fopen(filename,'wt');
    fprintf(fid,'%d\n\n',nSol);
    
    for k = 1:nSol
        for l = 1:lSol
            fprintf(fid,'%.16f\t',solution(l,k));
            fprintf(fid,'%.16f\n',0);
        end
        fprintf(fid,'\n');
    end
    fclose(fid); 
end