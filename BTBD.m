%%% Housekeeping
clc; clear; close all;
rng(41)

%%% To Do
% use fix bond react.
% update readme.

%%% Notation
% r - position vector
% r12 - vector pointing from position of bead 1 to bead 2
% o, os, oi - origami, list of origamis, origami index
% b, bs, bi - block, list of blocks, block index
% nvar - number of var
% var.n or n_var - number of beads in var
% ib - bead index within block
% io - bead index within entire origami
% iu - bead index within the universe
% connection - permenant (usually scaffold) bond
% linker - switchable (usually sticky end) bond


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Heart %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% read input
inFile = "./designs/triarm_ds3.txt";
[p,os,linker_types,pot_types,apot_types,nABAtype] = read_input(inFile);

%%% set output
outFold = "/Users/dduke/Files/block_tether/network/experiment/active/";
nsim = 1;

%%% create output folder 
mkdir(outFold)

%%% loop over simulations
for i = 1:nsim
    if nsim == 1
        simFold = outFold;
    else
        simFold = outFold + "sim" + ars.fstring(i,2,0,"R","zero") + "/";
        mkdir(simFold)
    end

    %%% initialize positions
    os = init_positions(os,p);

    %%% write lammps simulation geometry file
    geoFile = simFold + "geometry.in";
    geoVisFile = simFold + "geometry_vis.in";
    compose_geo(geoFile,geoVisFile,p.dbox,nABAtype,os);
    
    %%% write lammps input file
    inputFile = simFold + "lammps.in";
    write_input(inputFile,p,linker_types,pot_types,apot_types)
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% File Readers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% read parameters from input file
function p = read_parameters(inFile)

    %%% define parameters with their default values
    p_input = struct( ...
        'nstep',        NaN, ...        % steps         - number of production steps (required)
        'nstep_relax',  1e5, ...        % steps         - number of relaxation/shrink steps
        'dump_every',   1e4, ...        % steps         - number of steps between dumps
        'dbox',         NaN, ...        % nm            - periodic boundary diameter
        'shrink_ratio', 1,   ...        % none          - box compression (final/initial)
        'dt',           NaN, ...        % ns            - integration time step
        'verlet_skin',  4,   ...        % nm            - width of neighbor list skin
        'neigh_every',  1e1, ...        % steps         - how often to consider updating neighbor list
        'react_every',  1e1, ...        % steps         - how often to check for linker hybridization
        'T',            300, ...        % K             - temperature
        'sigma',        NaN, ...        % nm            - WCA distance parameter
        'epsilon',      NaN, ...        % kcal/mol      - WCA energy parameter
        'r12_bead',     NaN, ...        % nm            - block bead separation in xy-plane
        'r12_helix',    NaN, ...        % nm            - block bead separation along z-axis
        'U_strained',   10   ...        % kcal/mol      - max energy for initialized bonds and angles
    );

    %%% open file
    f = fopen(inFile, 'r');
    if f == -1
        error("Could not open file.");
    end

    %%% loop over lines
    while ~feof(f)
        line = strtrim(fgetl(f));

        %%% clean up line
        extract = split(extractBefore([line,'%'],'%'));
        extract = extract(~cellfun('isempty',extract));
        if length(extract) ~= 2
            continue
        end
   
        %%% set parameter
        if ismember(extract{1}, fieldnames(p_input))
            p_input.(extract{1}) = str2double(extract{2});
        end
    end
    fclose(f);

    %%% check for missing parameters
    missing = cellfun(@(x) isnumeric(x) && any(isnan(x)), struct2cell(p_input));
    if any(missing)
        fields = fieldnames(p_input);
        error('Missing required parameters: %s', strjoin(fields(missing), ', '));
    end

    %%% parameters
    p = parameters(p_input);
end


%%% read input file and create corresponding origami objects
function [p,os,linker_types,pot_types,apot_types,nABAtype] = read_input(inFile)

    %%% read parameters
    p = read_parameters(inFile);

    %%% open file
    f = fopen(inFile, 'r');
    if f == -1
        error("Could not open file.");
    end

    %%% initialize dictionaries
    pot_types = dictionary();
    apot_types = dictionary();
    block_templates = dictionary();
    origami_templates = dictionary();
    origami_counts = dictionary();
    linker_types = dictionary();

    %%% interpret on and off
    status = dictionary("on",1,"off",0);

    %%% loop over lines
    while ~feof(f)
        line = strtrim(fgetl(f));

        %%% discard empty lines and whole-line comments
        if isempty(line) || startsWith(line, '%')
            continue
        end

        %%% extract values before comment
        extract = split(extractBefore([line,'%'],'%'));
        extract = extract(~cellfun('isempty',extract));
        if length(extract) < 2
            continue
        end

        %%% log input
        switch extract{1}

            %%% define bond potential
            case 'bond_pot'
                style = string(extract{3});
                r12_eq = str2double(extract{4});
                params = str2double(extract(5:end));
                index = 1 + numEntries(pot_types);
                pot_types(extract{2}) = bond_pot(style,r12_eq,params,index);

            %%% define angle potential
            case 'angle_pot'
                theta_eq = str2double(extract{3});
                k_theta = str2double(extract{4});
                index = 1+numEntries(apot_types);
                apot_types(extract{2}) = angle_pot(theta_eq,k_theta,index);

            %%% define block
            case 'block'
                if strcmp(extract{3},'copy')
                    block_templates(extract{2}) = block_templates(extract{4});
                else
                    pattern_label = convertCharsToStrings(extract{3});
                    n_helix = str2double(extract{4});
                    block_templates(extract{2}) = block(pattern_label,n_helix,p);
                end

            %%% add patch to block
            case 'patch'
                name = extract{2};
                x = str2double(extract{4});
                y = str2double(extract{5});
                z = str2double(extract{6});
                block_templates(extract{3}) = block_templates(extract{3}).add_patch(name,x,y,z);

                %%% keyword value pairs
                line_index = 7;
                while true
                    if length(extract) < line_index
                        break
                    else
                        switch extract{line_index}
                            case 'status'
                                block_templates(extract{3}).status(end) = status(extract{line_index+1});
                                line_index = line_index + 2;
                            otherwise
                                error("Unknown patch keyword: " + extract{line_index})
                        end
                    end
                end

            %%% initialize origami
            case 'origami'
                index = 1 + numEntries(origami_templates);
                origami_templates(extract{2}) = origami(index);
                origami_counts(extract{2}) = str2double(extract{3});

                %%% keyword value pairs
                line_index = 4;
                while true
                    if length(extract) < line_index
                        break
                    else
                        switch extract{line_index}
                            case 'copy'
                                origami_templates(extract{2}) = origami_templates(extract{line_index+1});
                                line_index = line_index + 2;
                            case 'rigid'
                                origami_templates(extract{2}).rigid = string(extract{line_index+1});
                                line_index = line_index + 2;
                            otherwise
                                error("Unknown origami keyword: " + extract{line_index})
                        end
                    end
                end

            %%% initialize linker
            case 'linker'
                pot_index = pot_types(extract{3}).index;
                r12_cut = str2double(extract{4});
                index = 1 + numEntries(linker_types);
                linker_types(extract{2}) = linker(pot_index,r12_cut,index);
            
           %%% add features to origamis and linkers
            otherwise

                %%% add feature to origami
                if numEntries(origami_templates) > 0 && isKey(origami_templates,extract{1})
                    switch extract{2}

                        %%% add blocks
                        case 'blocks'
                            for ei = 3:length(extract)
                                b = block_templates(extract{ei});
                                origami_templates(extract{1}) = origami_templates(extract{1}).add_block(b);
                            end

                        %%% add connection
                        case 'conn'
                            bi1 = str2double(extract{3});
                            patch1 = extract{4};
                            bi2 = str2double(extract{5});
                            patch2 = extract{6};
                            pot = pot_types(extract{7});
                            origami_templates(extract{1}) = origami_templates(extract{1}).add_conn(bi1,patch1,bi2,patch2,pot);
                        
                            %%% keyword value pairs
                            line_index = 8;
                            while true
                                if length(extract) < line_index
                                    break
                                else
                                    switch extract{line_index}
                                        case 'status'
                                            origami_templates(extract{1}).conns_status(end) = status(extract{line_index+1});
                                            line_index = line_index + 2;
                                        otherwise
                                            error("Unknown origami connection keyword: " + extract{line_index})
                                    end
                                end
                            end

                        %%% add angle
                        case 'angle'
                            bi1 = str2double(extract{3});
                            patch11 = extract{4};
                            patch12 = extract{5};
                            bi2 = str2double(extract{6});
                            patch21 = extract{7};
                            patch22 = extract{8};
                            apot = apot_types(extract{9});
                            origami_templates(extract{1}) = origami_templates(extract{1}).add_angle(bi1,patch11,patch12,bi2,patch21,patch22,apot);

                            %%% keyword value pairs
                            line_index = 10;
                            while true
                                if length(extract) < line_index
                                    break
                                else
                                    switch extract{line_index}
                                        case 'theta_init'
                                            theta_init = str2double(extract{line_index+1});
                                            origami_templates(extract{1}).angles_theta_init(end) = theta_init;
                                            line_index = line_index + 2;
                                        case 'axis_init'
                                            axis_init = [str2double(extract{line_index+1});str2double(extract{line_index+2});str2double(extract{line_index+3})];
                                            origami_templates(extract{1}).angles_axis_init(:,end) = axis_init;
                                            line_index = line_index + 4;
                                        case 'status'
                                            origami_templates(extract{1}).angles_status(end) = status(extract{line_index+1});
                                            line_index = line_index + 2;
                                        otherwise
                                            error("Unknown origami angle keyword: " + extract{line_index})
                                    end
                                end
                            end

                        %%% error
                        otherwise
                            error("Unknown origami parameter: " + extract{2})
                    end

                %%% add feature to linker
                elseif numEntries(linker_types) > 0 && isKey(linker_types,extract{1})
                    switch extract{2}

                        %%% define 5' end
                        case '5p'
                            li = linker_types(extract{1}).index;
                            bi = extract{4};
                            patch = extract{5};
                            origami_templates(extract{3}) = origami_templates(extract{3}).add_linker(1,li,bi,patch);
                        
                        %%% define 3' end
                        case '3p'
                            li = linker_types(extract{1}).index;
                            bi = extract{4};
                            patch = extract{5};
                            origami_templates(extract{3}) = origami_templates(extract{3}).add_linker(0,li,bi,patch);

                        %%% add angle
                        case 'angle'
                            linker_types(extract{1}).apot_index = apot_types(extract{3}).index;

                            %%% keyword value pairs
                            line_index = 4;
                            while true
                                if length(extract) < line_index
                                    break
                                else
                                    switch extract{line_index}
                                        case 'theta_min'
                                            linker_types(extract{1}).theta_min = str2double(extract{line_index+1});
                                            line_index = line_index + 2;
                                        case 'theta_max'
                                            linker_types(extract{1}).theta_max = str2double(extract{line_index+1});
                                            line_index = line_index + 2;
                                        otherwise
                                            error("Unknown linker keyword: " + extract{line_index})
                                    end
                                end
                            end
                        
                        %%% error
                        otherwise
                            error("Unknown linker parameter: " + extract{2})
                    end

                %%% error
                elseif ~ismember(extract{1},fieldnames(p))
                    error("Unknown system parameter: " + extract{1})
                end
        end
    end
    fclose(f);

    %%% create origamis
    os = origami.empty;
    for o_name = keys(origami_templates)'
        count = origami_counts(o_name);
        os(length(os)+1:length(os)+count) = origami_templates(o_name);
    end

    %%% count atom, bond, angle types
    natomType = 2 + numEntries(linker_types)*2;
    nbondType = 1 + numEntries(pot_types);
    nangleType = numEntries(apot_types);
    nABAtype = [natomType,nbondType,nangleType];
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% File Writers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% write lammps geometry file
function compose_geo(geoFile,geoVisFile,dbox,nABAtype,os)

    %%% initialize mass info
    masses = ones(1,nABAtype(1));
    mass_patch = 0.01;
    masses(2) = mass_patch;

    %%% count origami types
    norigami_type = length(unique([os.index]));

    %%% initialize universal index map
    get_iu = zeros(length(os),max([os.n]));

    %%% compile atom info
    atoms = zeros(5,0);
    atoms_vis = zeros(5,0);
    atom_count = 0;
    rigid_count = 0;
    for oi = 1:length(os)
        for io = 1:os(oi).n

            %%% determine bead status
            bi = os(oi).get_bi(io);
            ib = os(oi).get_ib(io);
            if os(oi).bs(bi).status(ib) == 1
                atom_count = atom_count + 1;
                get_iu(oi,io) = atom_count;

                %%% increment rigid body count if new body
                if os(oi).rigid == "origami"
                    if io == 1 
                        rigid_count = rigid_count + 1;
                    end
                elseif os(oi).rigid == "block"
                    if os(oi).get_ib(io) == 1
                        rigid_count = rigid_count + 1;
                    end
                else
                    error("Unknown rigid type")
                end

                %%% atom info
                atoms(1,atom_count) = rigid_count;
                if os(oi).is_patch(io)
                    atoms(2,atom_count) = 2;
                else
                    atoms(2,atom_count) = 1;
                end
                atoms(3:5,atom_count) = os(oi).r(:,io);

                %%% visualization info
                atoms_vis(1,atom_count) = oi;
                if os(oi).is_patch(io)
                    atoms_vis(2,atom_count) = norigami_type + 1;
                else
                    atoms_vis(2,atom_count) = os(oi).index;
                end
                atoms_vis(3:5,atom_count) = os(oi).r(:,io);
            end
        end

        %%% edit types for linkers
        for lio = 1:os(oi).nlink5
            io = os(oi).link5s_io(lio);
            iu = get_iu(oi,io);
            if iu == 0
                error("Linker bead does not exist.")
            end
            li = os(oi).link5s_index(lio);
            atoms(2,iu) = 1 + 2*li;
            atoms_vis(2,iu) = norigami_type + 1 + li;
            if os(oi).is_patch(io)
                masses(1+2*li) = mass_patch;
            end
        end
        for lio = 1:os(oi).nlink3
            io = os(oi).link3s_io(lio);
            iu = get_iu(oi,io);
            if iu == 0
                error("Linker bead does not exist.")
            end
            li = os(oi).link3s_index(lio);
            atoms(2,iu) = 2 + 2*li;
            atoms_vis(2,iu) = norigami_type + 1 + li;
            if os(oi).is_patch(io)
                masses(2+2*li) = mass_patch;
            end
        end
    end

    %%% connection bonds
    bonds = zeros(3,0);
    bond_count = 0;
    for oi = 1:length(os)
        for ci = 1:length(os(oi).conns_pot)
            if os(oi).conns_status(ci) == 1
                bond_count = bond_count + 1;
                type = os(oi).conns_pot(ci).index;
                iu1 = get_iu( oi, os(oi).get_io( os(oi).conns_bis(1,ci), os(oi).conns_ibs(1,ci) ) );
                iu2 = get_iu( oi, os(oi).get_io( os(oi).conns_bis(2,ci), os(oi).conns_ibs(2,ci) ) );
                if iu1 == 0 || iu2 == 0
                    error("Connection bead does not exist.")
                end
                bonds(1,bond_count) = type;
                bonds(2,bond_count) = iu1;
                bonds(3,bond_count) = iu2;
            end
        end
    end

    %%% rigid bonds for angles
    for oi = 1:length(os)
        for ai = 1:length(os(oi).angles_apot)
            if os(oi).angles_status(ai) == 1
                iu1 = get_iu( oi, os(oi).get_io( os(oi).angles_bis(1,ai), os(oi).angles_ibs(1,ai) ) );
                iu2 = get_iu( oi, os(oi).get_io( os(oi).angles_bis(2,ai), os(oi).angles_ibs(2,ai) ) );
                if iu1 == 0 || iu2 == 0
                    error("Angle bead does not exist.")
                end
                bond_count = bond_count + 1;
                bonds(1,bond_count) = nABAtype(2);
                bonds(2,bond_count) = iu1;
                bonds(3,bond_count) = iu2;
                iu3 = get_iu( oi, os(oi).get_io( os(oi).angles_bis(3,ai), os(oi).angles_ibs(3,ai) ) );
                iu4 = get_iu( oi, os(oi).get_io( os(oi).angles_bis(4,ai), os(oi).angles_ibs(4,ai) ) );
                if iu3 == 0 || iu4 == 0
                    error("Angle bead does not exist.")
                end
                bond_count = bond_count + 1;
                bonds(1,bond_count) = nABAtype(2);
                bonds(2,bond_count) = iu3;
                bonds(3,bond_count) = iu4;
            end
        end
    end

    %%% compile angle info
    angles = zeros(4,0);
    angle_count = 0;
    for oi = 1:length(os)
        for ai = 1:length(os(oi).angles_apot)
            if os(oi).angles_status(ai) == 1
                type = os(oi).angles_apot(ai).index;
                iu1 = get_iu( oi, os(oi).get_io( os(oi).angles_bis(1,ai), os(oi).angles_ibs(1,ai) ) );
                iu2 = get_iu( oi, os(oi).get_io( os(oi).angles_bis(2,ai), os(oi).angles_ibs(2,ai) ) );
                iu3 = get_iu( oi, os(oi).get_io( os(oi).angles_bis(3,ai), os(oi).angles_ibs(3,ai) ) );
                iu4 = get_iu( oi, os(oi).get_io( os(oi).angles_bis(4,ai), os(oi).angles_ibs(4,ai) ) );
                angle_count = angle_count + 1;
                angles(1,angle_count) = type;
                angles(2,angle_count) = iu1;
                angles(3,angle_count) = iu2;
                angles(4,angle_count) = iu3;
                angle_count = angle_count + 1;
                angles(1,angle_count) = type;
                angles(2,angle_count) = iu2;
                angles(3,angle_count) = iu3;
                angles(4,angle_count) = iu4;
            end
        end
    end

    %%% write simulation geometry file
    ars.writeGeo(geoFile,dbox,atoms,bonds,angles,masses=masses,nbondType=nABAtype(2),nangleType=nABAtype(3));

    %%% write visualization geometry file
    ars.writeGeo(geoVisFile,dbox,atoms_vis,bonds,angles)
end


%%% write lammps input file
function write_input(inputFile,p,linker_types,potential_types,angle_types)

    %%% get counts
    nlinker = numEntries(linker_types);
    natomType = 2 + nlinker*2;
    nbondType = 1 + numEntries(potential_types);
    nangleType = numEntries(angle_types);

    %%% calculate communication cutoff
    U_max = 12;
    comm_cutoff = p.r12_cut_WCA;
    if numEntries(potential_types) > 0
        for p_name = keys(potential_types)'
            r12_max = potential_types(p_name).calc_separation(U_max);
            comm_cutoff = max([comm_cutoff,r12_max]);
        end
    end
    if numEntries(linker_types) > 0
        for l_name = keys(linker_types)'
            r12_cut = linker_types(l_name).r12_cut;
            comm_cutoff = max([comm_cutoff,r12_cut]);
        end
    end

    %%% open file
    f = fopen(inputFile,'w');
    
    %%% header
    fprintf(f,strcat(...
        "\n#------ Begin Input ------#\n",...
        "# Written by BTBD.m\n\n"));

    %%% basic setup
    fprintf(f,strcat(...
        "## Geometry\n",...
        "units           nano\n",...
        "dimension       3\n",...
        "boundary        p p p\n",...
        "atom_style      molecular\n",...
        "read_data       geometry.in &\n",...
        "                extra/bond/per/atom 1 &\n",...
        "                extra/special/per/atom 1\n\n"));

    %%% neighbor list
    fprintf(f,strcat(...
        "## Parameters\n",...
        "neighbor        ", ars.fstring(p.verlet_skin,0,2), " bin\n",...
        "neigh_modify    every ", num2str(p.neigh_every), "\n",...
        "neigh_modify    exclude molecule/intra all\n"));

    %%% pairwise interactions
    fprintf(f,strcat(...
        "pair_style      hybrid/overlay lj/cut ", ars.fstring(p.r12_cut_WCA,0,2), " zero ", ars.fstring(comm_cutoff,0,2), "\n",...
        "pair_coeff      * * zero\n",...
        "pair_coeff      1 1 lj/cut ", ars.fstring(p.epsilon,0,2), " ", ars.fstring(p.sigma,0,2), "\n",...
        "pair_modify     pair lj/cut shift yes\n"));

    %%% bonds
    if nbondType == 1
        fprintf(f,strcat(...
            "bond_style      zero\n",...
            "bond_coeff      1\n"));
    else
        fprintf(f,...
            "bond_style      hybrid");
        styles = unique(["zero",values(potential_types).style]);
        for si = 1:length(styles)
            fprintf(f,strcat(" ", styles{si}));
        end
        fprintf(f,"\n");
        for p_name = keys(potential_types)'
            potential_types(p_name).write_potential(f);
        end
        fprintf(f,strcat(...
            "bond_coeff      ", num2str(nbondType), " zero\n"));
    end

    %%% angles
    if nangleType > 0
        fprintf(f,...
            "angle_style     harmonic\n");
        for a_name = keys(angle_types)'
            fprintf(f,strcat(...
                "angle_coeff     ", num2str(angle_types(a_name).index), " ", ars.fstring(angle_types(a_name).k_theta/2,0,2), " ", ars.fstring(angle_types(a_name).theta_eq,0,2), "\n"));
        end
    end

    %%% linker group
    if nlinker > 0
        fprintf(f,strcat(...
            "group           linker type 3:", num2str(natomType), "\n"));
    end
    fprintf(f,"\n");

    %%% thermostat
    fprintf(f,strcat(...
        "## Thermostat\n",...
        "fix             tstat all rigid/nve molecule langevin ", num2str(p.T), " ", num2str(p.T), " ", num2str(p.dt*10), " 37\n",...
        "thermo          ", num2str(p.dump_every), "\n\n"));

    %%% relaxation
    if p.nstep_relax > 0
        fprintf(f,strcat(...
            "## Relaxation\n",...
            "timestep        ", num2str(p.dt/10), "\n"));
        if p.shrink_ratio ~= 1
            fprintf(f,strcat(...
                "fix             shrink all deform 1 &\n",...
                "                x final ", ars.fstring(-p.dbox/2*p.shrink_ratio,0,2), " ", ars.fstring(p.dbox/2*p.shrink_ratio,0,2), " &\n",...
                "                y final ", ars.fstring(-p.dbox/2*p.shrink_ratio,0,2), " ", ars.fstring(p.dbox/2*p.shrink_ratio,0,2), " &\n",...
                "                z final ", ars.fstring(-p.dbox/2*p.shrink_ratio,0,2), " ", ars.fstring(p.dbox/2*p.shrink_ratio,0,2), "\n"));
        end
        fprintf(f,strcat(...
            "run             ", num2str(p.nstep_relax), "\n"));
        if p.shrink_ratio ~= 1
            fprintf(f,strcat(...
                "unfix           shrink\n"));
        end
        fprintf(f,strcat(...
            "reset_timestep  0\n\n"));
    end

    %%% reactions
    if nlinker > 0
        fprintf(f,"## Reactions\n");
        for l_name = keys(linker_types)'
            li = linker_types(l_name).index;
            fix_name = strcat("linker",num2str(li));
            r12_cut = linker_types(l_name).r12_cut;
            bond_type = linker_types(l_name).pot_index;
            angle_type = linker_types(l_name).apot_index;
            theta_min = linker_types(l_name).theta_min;
            theta_max = linker_types(l_name).theta_max;
            write_bond_create(f,fix_name,"linker",p.react_every,1+li*2,2+li*2,r12_cut,bond_type,angle_type,theta_min,theta_max);
        end
        fprintf(f,"\n");
    end
    
    %%% updates
    if p.nstep > 0
        fprintf(f,strcat(...
            "## Updates\n",...
            "dump            dumpT all custom ", num2str(p.dump_every), " trajectory.dat id mol xs ys zs\n",...
            "dump_modify     dumpT sort id\n",...
            "compute         compB1 all bond/local dist engpot\n",...
            "compute         compB2 all property/local btype batom1 batom2\n",...
            "dump            dumpB all local ", num2str(p.dump_every), " dump_bonds.txt index c_compB1[1] c_compB1[2] c_compB2[1] c_compB2[2] c_compB2[3]\n"));
        if nangleType > 0
            fprintf(f,strcat(...
			"compute         compA1 all angle/local theta eng\n",...
			"compute         compA2 all property/local atype aatom1 aatom2 aatom3\n",...
		    "dump            dumpA all local ", num2str(p.dump_every), " dump_angles.dat index c_compA1[1] c_compA1[2] c_compA2[1] c_compA2[2] c_compA2[3] c_compA2[4]\n"));
        end
        fprintf(f,"\n");
    end

    %%% production
    if p.nstep > 0
        fprintf(f,strcat(...
            "## Production\n",...
            "timestep        ", num2str(p.dt), "\n",...
            "run             ", num2str(p.nstep), "\n"));
    end
    
    %%% finalize and close file
    fprintf(f,"\n#------- End Input -------#\n\n");
    fclose(f);
end


%%% write fix bond create command
function write_bond_create(f,fix_name,group_name,react_every,atomType1,atomType2,r12_cut,bond_type,angle_type,theta_min,theta_max)
    fprintf(f,strcat(...
        "fix             ", fix_name, " ", group_name, " bond/create"));
    if angle_type ~= 0
        fprintf(f,"/angle");
    end
    fprintf(f,strcat(...
        " ", num2str(react_every), " ",...
        num2str(atomType1), " ",num2str(atomType2), " ",...
        ars.fstring(r12_cut,0,2), " ", num2str(bond_type), " "));
    if angle_type ~= 0 
        fprintf(f,strcat(...
            "atype ", num2str(angle_type), " ",...
            "aconstrain ", num2str(theta_min), " ", num2str(theta_max), " "));
    end
    fprintf(f,strcat(...
        "iparam 1 ", num2str(atomType1), " ",...
        "jparam 1 ", num2str(atomType2), "\n"));
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% initilize positions of the entire system
function os = init_positions(os,p)
    max_attempts = 10;

    %%% system attempt loop
    attempts = 0;
    while true

        %%% initialize avoided positions
        r_other = [];

        %%% check system attempts
        if attempts == max_attempts
            error("Could not place origamis.")
        end

        %%% loop over origamis
        for oi = 1:length(os)

            %%% place origami
            disp(strcat("Initializing origami ",num2str(oi),"..."))
            [os(oi),failed,r_other] = os(oi).init_positions(p,r_other);

            %%% stop looping through origamis if one failed
            if failed == true
                break
            end
        end

        %%% check if origami loop succeeded
        if failed == true
            attempts = attempts + 1;
            continue
        end

        %%% system successfully initiated
        fprintf("Initialization complete.\n")
        break
    end
end

