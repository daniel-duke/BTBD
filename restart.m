%%% read old lammps file and write new one to continue the simulation
function restart(nstep,replace,nsim,outFold)

    %%% set arguments
    arguments
        nstep double
        replace double = 0
        nsim double = 1
        outFold char = pwd
    end

    %%% loop over simulations
    for i = 1:nsim

        %%% single simulation
        if nsim == 1
            simFold = outFold;

        %%% multiple simulations
        else
            if i == 1; ndigit_nsim = floor(log10(nsim))+1; end
            if i > 1; p.rseed_lmp = p.rseed_lmp + 1; end
            simFileName = "sim" + ars.fstring(i,ndigit_nsim,0,"R","zero");
            simFold = fullfile(outFold,simFileName);
        end

        %%% get current run number
        lammpsFile = fullfile(simFold,"lammps.in");
        run = get_run_number(lammpsFile);

        %%% save copy of old file
        if ~replace
            oldLammpsFileName = "lammps_run" + run + ".in";
            copyfile(lammpsFile,fullfile(simFold,oldLammpsFileName))
            run = run + 1;
        end

        %%% rewrite lammps file
        edit_lammps(lammpsFile,nstep,run)
    end
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% File Handlers %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% extract run number from lammps file
function run = get_run_number(lammpsFile)

    %%% open file
    ars.checkFileExist(lammpsFile,"lammps");
    f = fopen(lammpsFile,'r');

    %%% read file
    content = textscan(f, '%s', 'Delimiter', '\n', 'Whitespace', '');
    content = content{1};
    fclose(f);

    %%% find run number
    for i = 1:length(content)
        if startsWith(content{i}, '# Run')
            run = str2double(content{i}(7:end));
            return
        end
    end
end


%%% adjust lammps script for restarting the simulation
function edit_lammps(lammpsFile,nstep,run)

    % read old lammps file
    f = fopen(lammpsFile,'r');
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

        %%% adjust run number
        if startsWith(content_in{i_in}, '# Run')
            content_out{i_out} = "# Run " + num2str(run);

        %%% set how to read geometry
        elseif startsWith(content_in{i_in}, '## Geometry')
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
    f = fopen(lammpsFile,'w');
    fprintf(f, '%s\n', content_out{:});
    fclose(f);
end
