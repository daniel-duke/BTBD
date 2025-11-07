%%% Housekeeping
clc; clear; close all;

%%% To Do
% write the script


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Heart %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% parameters
nstep = 1E6;

%%% set output
outFold = "/Users/dduke/Data/triarm/experiment/active/";
nsim = 1;

%%% loop over simulations
for i = 1:nsim
    if nsim == 1
        simFold = outFold;
    else
        simFold = outFold + "sim" + ars.fstring(i,2,0,"R","zero") + "/";
    end

    %%% rewrite lammps input file
    inputFile = simFold + "lammps.in";
    edit_input(inputFile, nstep)
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% File Handlers %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% read old lammps file and write edited version
function edit_input(inFile, nstep)

    %%% check file exists
    if ~isfile(inFile)
        error("Could not find input file: " + inFile);
    end

    % read old lammps file
    f = fopen(inFile, 'r');
    content_in = textscan(f, '%s', 'Delimiter', '\n', 'Whitespace', '');
    content_in = content_in{1};
    fclose(f);

    % parse the content and write edited content
    content_out = cell(1);
    i_in = 1;
    i_out = 1;

    %%% parse content
    while i_in <= numel(content_in)
        content_out{i_out} = content_in{i_in};

        %%% set how to read geometry
        if startsWith(content_in{i_in}, '## Geometry')
            content_out{i_out+1} = 'read_restart    restart_binary2.out';
            i_out = i_out + 1;

            %%% skip the appropriate number of lines
            if startsWith(content_in{i_in+1}, 'read_restart')
                i_in = i_in + 1;
            else
                i_in = i_in + 5;
            end
        
        %%% remove relaxation
        elseif startsWith(content_in{i_in}, '## Relaxation')
            content_out{i_out} = '';
            i_out = i_out - 1;
            i_in = i_in + 4;
        
        %%% run a single step before dumps
        elseif startsWith(content_in{i_in}, '## Production')
            if ~startsWith(content_in{i_in+2}, 'run')
                content_out{i_out+1} = content_in{i_in+1};
                content_out{i_out+2} = 'run             1';
                i_out = i_out + 2;
                i_in = i_in + 1;
            end

        %%% adjust run time
        elseif startsWith(content_in{i_in}, '## Go Time')
            line = "run             " + num2str(nstep-1);
            content_out{i_out+1} = convertStringsToChars(line);
            i_out = i_out + 1;
            i_in = i_in + 1;
        end

        %%% next line
        i_in = i_in + 1;
        i_out = i_out + 1;
    end

    %%% write edited lammps file
    f = fopen(inFile, 'w');
    fprintf(f, '%s\n', content_out{:});
    fclose(f);
end
