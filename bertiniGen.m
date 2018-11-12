%===================== Generate Bertini input files =======================

clear; close all; clc;
addpath(genpath('src'));
cur_path=pwd;

%% Load and define parameters
load('Data/parameters_sep.mat');
rho = p.rho;
drho = .0100000;

%% Input interface
inputMsg = ['Which kind of input file do you want?\n'...
              '1 -- Regular input;  '...
              '2 -- Tracking with rho;  '...  
              '3 -- Homotopy tracking;\n'...
              '4 -- With linearized equations;  '...
              '5 -- With parameter boundaries\n'...
              'Your selection: \n'];
getHom = input(inputMsg);

%% Write input files
% multiple grid files
inputMsg = ['Do you want input files for multiple grids?\n'...
              '1 -- Yes;  '... 
              'Else -- No.\n'...
              'Your selection: \n'];
multi_inputs = input(inputMsg);
if multi_inputs == 1
    for k = 0:5

        switch getHom
            case 1
                N = 5*2^k;
                MyPath = sprintf('bertiniFile/Sep/input_%d', N);
                MyFilename=sprintf('%s/GS_input_%d', MyPath, N);
                userhom = 0;

            case 2
                N = 5*2^k;
                MyPath = sprintf('bertiniFile/input_rho_%d', N);
                MyFilename=sprintf('%s/GS_input_track_rho_%d', MyPath, N);
                userhom = 1;

            case 3
                N = 5*2^k;
                MyPath = sprintf('bertiniFile/input_hom_%d', N);
                MyFilename=sprintf('%s/GS_input_hom_%d', MyPath, N);
                userhom = 1;
                
            case 4
                N = 5*2^k;
                MyPath = sprintf('bertiniFile/input_lin_%d', N);
                MyFilename=sprintf('%s/GS_input_lin_%d', MyPath, N);
                userhom = 0;
                
            case 5
                N = 5*2^k;
                MyPath = sprintf('bertiniFile/input_bdry_%d', N);
                MyFilename=sprintf('%s/GS_input_bdry_%d', MyPath, N);
                userhom = 0;

            otherwise
            error('Invalid input file name!')
        end

        if ~isdir(MyPath)
            mkdir(MyPath);
        end

        switch getHom
            case 1
                GS_input(N,MyFilename,userhom,p);

            case 2
                GS_input_track_rho(N,MyFilename,userhom,p);

            case 3
                GS_input_hom(N,MyFilename,userhom,p);

            case 4
                GS_input_lin(N,MyFilename,userhom,p);
                
            case 5
                GS_input_bdry(N,MyFilename,userhom,p);
        end
    end
end

% track rho
inputMsg = ['Do you want to track rho for N=5*2^k?\n'...
              '1 (i.e. 10-160) -- Yes;  '... 
              '0 -- No tracking.\n'...
              'Your selection: \n'];
track_rho = input(inputMsg);
if track_rho
    %parpool;
    for k = 1:5
        MyPath = sprintf('bertiniFile/input_rho_%d',5*2^k);
        MyFilename=sprintf('%s/GS_input_rho_%d', MyPath,5*2^k);
        
        if ~isdir(MyPath)
            mkdir(MyPath);
        end
        delete(MyFilename);
        GS_input_track_rho(5*2^k,MyFilename,1,p);
        
        % some terminal commands
%         str = sprintf('cp start_gen %s/start', MyPath);
%         system(str);
%         str = sprintf('cd %s/%s', cur_path, MyPath);
%         system(str);
%         str = sprintf('nohup /usr/local/bin/bertini %s', MyFilename);
%         system(str);
    end
end

