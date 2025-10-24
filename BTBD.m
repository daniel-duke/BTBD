%%% Housekeeping
clc; clear; close all;
rng(42)

%%% To Do
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
inFile = "./designs/triarm_ds3_se1.txt";
[p,os,ls,rs,pots,apots,dpots,nABADtype] = read_input(inFile);

%%% set output
outFold = "/Users/dduke/Files/triarm/experiment/active3/";
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

    %%% write lammps input file
    inputFile = simFold + "lammps.in";
    reactFold = simFold + "react/";
    write_input(inputFile,reactFold,p,nABADtype,ls,rs,pots,apots,dpots)

    %%% write lammps simulation geometry file
    geoFile = simFold + "geometry.in";
    geoVisFile = simFold + "geometry_vis.in";
    compose_geo(geoFile,geoVisFile,p.dbox,nABADtype,os,ls,rs);
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
        'U_strained',   10,  ...        % kcal/mol      - max energy for initialized bonds and angles
        'comm_cutoff',  0    ...        % nm            - communication_cutoff
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
function [p,os,ls,rs,pots,apots,dpots,nABADtype] = read_input(inFile)

    %%% read parameters
    p = read_parameters(inFile);

    %%% open file
    f = fopen(inFile, 'r');
    if f == -1
        error("Could not open file.");
    end

    %%% initialize dictionaries
    pots = dictionary();
    apots = dictionary();
    dpots = dictionary();
    block_templates = dictionary();
    origami_templates = dictionary();
    origami_counts = dictionary();
    linkers = dictionary();
    reactions = dictionary();

    %%% define zero potential
    pots("zero") = bond_pot("zero",0,zeros(0,1),1);

    %%% keep track of atom types
    ti_count = 2;

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

        %%% read input
        switch extract{1}

            %%% define bond potential
            case 'bond_pot'
                style = string(extract{3});
                r12_eq = str2double(extract{4});
                params = str2double(extract(5:end));
                index = 1 + numEntries(pots);
                pots(extract{2}) = bond_pot(style,r12_eq,params,index);

            %%% define angle potential
            case 'angle_pot'
                style = string(extract{3});
                theta_eq = str2double(extract{4});
                params = str2double(extract(5:end));
                index = 1+numEntries(apots);
                apots(extract{2}) = angle_pot(style,theta_eq,params,index);

            %%% define dihedral potential
            case 'dihedral_pot'
                style = string(extract{3});
                theta_eq = str2double(extract{4});
                params = str2double(extract(5:end));
                index = 1+numEntries(dpots);
                dpots(extract{2}) = dihedral_pot(style,theta_eq,params,index);

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
                if length(extract) > 2 && strcmp(extract{3},'copy')
                    origami_templates(extract{2}) = origami_templates(extract{4});
                    origami_templates(extract{2}).label = string(extract{2});
                    line_index = 5;
                else
                    origami_templates(extract{2}) = origami(string(extract{2}));
                    line_index = 3;
                end

                %%% default origami count
                origami_counts(extract{2}) = 1;

                %%% keyword value pairs
                while true
                    if length(extract) < line_index
                        break
                    else
                        switch extract{line_index}
                            case 'count'
                                origami_counts(extract{2}) = str2double(extract{line_index+1});
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
                label = string(extract{2});
                pot_index = pots(extract{3}).index;
                r12_cut = str2double(extract{4});
                ti_start = ti_count + 1;
                linkers(extract{2}) = linker(label,pot_index,r12_cut,ti_start);
                ti_count = ti_count + 2;
            
            %%% initialize reaction
            case 'reaction'
                style = string(extract{3});
                params = extract(4:end);
                ti_start = ti_count + 1;
                reactions(extract{2}) = reaction(extract{2},style,params,ti_start);
                ti_count = ti_count + reactions(extract{2}).nti;

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
                            pot = pots(extract{7});
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
                            apot = apots(extract{9});
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

                        %%% add dihedral
                        case 'dihedral'
                            bi1 = str2double(extract{3});
                            patch1 = extract{4};
                            bi2 = str2double(extract{5});
                            patch2 = extract{6};
                            bi3 = str2double(extract{7});
                            patch3 = extract{8};
                            bi4 = str2double(extract{9});
                            patch4 = extract{10};
                            dpot = dpots(extract{11});
                            origami_templates(extract{1}) = origami_templates(extract{1}).add_dihedral(bi1,patch1,bi2,patch2,bi3,patch3,bi4,patch4,dpot);

                        %%% error
                        otherwise
                            error("Unknown origami parameter: " + extract{2})
                    end

                %%% add feature to linker
                elseif numEntries(linkers) > 0 && isKey(linkers,extract{1})
                    switch extract{2}

                        %%% add 5' end
                        case '5p'
                            o = origami_templates(extract{3});
                            ois = get_ois(origami_counts,extract{3});
                            bi = extract{4};
                            patch = string(extract(5));
                            linkers(extract{1}) = linkers(extract{1}).add_site(1,o,ois,bi,patch);

                        %%% add 3' end
                        case '3p'
                            o = origami_templates(extract{3});
                            ois = get_ois(origami_counts,extract{3});
                            bi = extract{4};
                            patch = string(extract(5));
                            linkers(extract{1}) = linkers(extract{1}).add_site(0,o,ois,bi,patch);

                        %%% add angle
                        case 'angle'
                            linkers(extract{1}).apot_index = apots(extract{3}).index;

                            %%% keyword value pairs
                            line_index = 4;
                            while true
                                if length(extract) < line_index
                                    break
                                else
                                    switch extract{line_index}
                                        case 'theta_min'
                                            linkers(extract{1}).theta_min = str2double(extract{line_index+1});
                                            line_index = line_index + 2;
                                        case 'theta_max'
                                            linkers(extract{1}).theta_max = str2double(extract{line_index+1});
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

                %%% add feature to reaction
                elseif numEntries(reactions) > 0 && isKey(reactions,extract{1})
                    switch extract{2}

                        %%% add 5' end
                        case '5p'
                            o = origami_templates(extract{3});
                            ois = get_ois(origami_counts,extract{3});
                            bi = extract{4};
                            patches = string(extract(5:end));
                            reactions(extract{1}) = reactions(extract{1}).add_site(1,o,ois,bi,patches);
                        
                        %%% add 3' end
                        case '3p'
                            o = origami_templates(extract{3});
                            ois = get_ois(origami_counts,extract{3});
                            bi = extract{4};
                            patches = string(extract(5:end));
                            reactions(extract{1}) = reactions(extract{1}).add_site(0,o,ois,bi,patches);
                    end

                %%% error
                elseif ~ismember(extract{1},fieldnames(p))
                    error("Unknown system parameter: " + extract{1})
                end
        end
    end
    fclose(f);

    %%% create origami array
    os = origami.empty;
    for o_name = keys(origami_templates)'
        ois = get_ois(origami_counts,o_name);
        os(ois) = origami_templates(o_name);
    end

    %%% create linker array
    ls = linker.empty;
    if numEntries(linkers) > 0
        ls = values(linkers)';
    end

    %%% create reaction array
    rs = reaction.empty;
    if numEntries(reactions) > 0
        rs = values(reactions)';
    end

    %%% count atom, bond, angle types
    natomType = ti_count;
    nbondType = numEntries(pots);
    nangleType = numEntries(apots);
    ndihedralType = numEntries(dpots);
    nABADtype = [natomType,nbondType,nangleType,ndihedralType];
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% File Writers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% write lammps geometry file
function compose_geo(geoFile,geoVisFile,dbox,nABADtype,os,ls,rs)
    disp("Writing geometry file...")

    %%% initialize mass info
    mass_patch = 0.01;
    masses = ones(1,nABADtype(1));
    masses(2:end) = mass_patch;

    %%% number origami types
    origami_types = dictionary();
    labels = unique([os.label]);
    norigamiType = length(labels);
    for i = 1:norigamiType
        origami_types(labels(i)) = i;
    end

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

                %%% molecule
                atoms(1,atom_count) = rigid_count;
                atoms_vis(1,atom_count) = oi;

                %%% type
                if ib <= os(oi).bs(bi).n_real
                    atoms(2,atom_count) = 1;
                    atoms_vis(2,atom_count) = origami_types(os(oi).label);
                else
                    atoms(2,atom_count) = 2;
                    atoms_vis(2,atom_count) = norigamiType + 1;
                end

                %%% position
                atoms(3:5,atom_count) = os(oi).r(:,io);
                atoms_vis(3:5,atom_count) = os(oi).r(:,io);
            end
        end
    end

    %%% edit types for linkers
    for li = 1:length(ls)
        for si = 1:ls(li).nsite5
            oi = ls(li).site5s_oi(si);
            bi = ls(li).site5s_bi(si);
            ib = ls(li).site5s_ib(si);
            io = os(oi).get_io(bi,ib);
            iu = get_iu(oi,io);
            ti = ls(li).ti_start;
            atoms(2,iu) = ti;
            atoms_vis(2,iu) = norigamiType + 1 + floor((ti-1)/2);
        end
        for si = 1:ls(li).nsite3
            oi = ls(li).site3s_oi(si);
            bi = ls(li).site3s_bi(si);
            ib = ls(li).site3s_ib(si);
            io = os(oi).get_io(bi,ib);
            iu = get_iu(oi,io);
            ti = ls(li).ti_start + 1;
            atoms(2,iu) = ti;
            atoms_vis(2,iu) = norigamiType + 1 + floor((ti-1)/2);
        end
    end

    %%% edit types for reactions
    for ri = 1:length(rs)
        for si = 1:rs(ri).nsite5
            oi = rs(ri).site5s_oi(si);
            bi = rs(ri).site5s_bi(si);
            for is = 1:rs(ri).n_site5
                ib = rs(ri).site5s_ibs(is,si);
                io = os(oi).get_io(bi,ib);
                iu = get_iu(oi,io);
                ti = rs(ri).tis_site5(is);
                atoms(2,iu) = ti;
                atoms_vis(2,iu) = norigamiType + 1 + floor((ti-1)/2);
            end
        end
        for si = 1:rs(ri).nsite3
            oi = rs(ri).site3s_oi(si);
            bi = rs(ri).site3s_bi(si);
            for is = 1:rs(ri).n_site3
                ib = rs(ri).site3s_ibs(is,si);
                io = os(oi).get_io(bi,ib);
                iu = get_iu(oi,io);
                ti = rs(ri).tis_site3(is);
                atoms(2,iu) = ti;
                atoms_vis(2,iu) = norigamiType + 1 + floor((ti-1)/2);
            end
        end
    end

    %%% connection bonds
    bonds = zeros(3,0);
    bond_count = 0;
    for oi = 1:length(os)
        for ci = 1:length(os(oi).conns_pot)
            if os(oi).conns_status(ci) == 1
                bond_type = os(oi).conns_pot(ci).index;
                iu1 = get_iu( oi, os(oi).get_io( os(oi).conns_bis(1,ci), os(oi).conns_ibs(1,ci) ) );
                iu2 = get_iu( oi, os(oi).get_io( os(oi).conns_bis(2,ci), os(oi).conns_ibs(2,ci) ) );
                if iu1 == 0 || iu2 == 0
                    error("Connection bead does not exist.")
                end
                bond_count = bond_count + 1;
                bonds(1,bond_count) = bond_type;
                bonds(2,bond_count) = iu1;
                bonds(3,bond_count) = iu2;
            end
        end
    end

    %%% rigid bonds for connection angles
    bond = zeros(3,1);
    bond_type = 1;
    for oi = 1:length(os)
        for ai = 1:length(os(oi).angles_apot)
            if os(oi).angles_status(ai) == 1
                iu1 = get_iu( oi, os(oi).get_io( os(oi).angles_bis(1,ai), os(oi).angles_ibs(1,ai) ) );
                iu2 = get_iu( oi, os(oi).get_io( os(oi).angles_bis(2,ai), os(oi).angles_ibs(2,ai) ) );
                if iu1 == 0 || iu2 == 0
                    error("Angle bead does not exist.")
                end
                bond(1) = bond_type;
                bond(2) = iu1;
                bond(3) = iu2;
                if ~any(all(bonds == bond, 1))
                    bond_count = bond_count + 1;
                    bonds(:,bond_count) = bond;
                end
                iu3 = get_iu( oi, os(oi).get_io( os(oi).angles_bis(3,ai), os(oi).angles_ibs(3,ai) ) );
                iu4 = get_iu( oi, os(oi).get_io( os(oi).angles_bis(4,ai), os(oi).angles_ibs(4,ai) ) );
                if iu3 == 0 || iu4 == 0
                    error("Angle bead does not exist.")
                end
                bond(1) = bond_type;
                bond(2) = iu3;
                bond(3) = iu4;
                if ~any(all(bonds == bond, 1))
                    bond_count = bond_count + 1;
                    bonds(:,bond_count) = bond;
                end
            end
        end
    end

    %%% rigid bonds for connection dihedrals
    bond = zeros(3,1);
    bond_type = 1;
    for oi = 1:length(os)
        for di = 1:length(os(oi).dihedrals_dpot)
            iu1 = get_iu( oi, os(oi).get_io( os(oi).dihedrals_bis(1,di), os(oi).dihedrals_ibs(1,di) ) );
            iu2 = get_iu( oi, os(oi).get_io( os(oi).dihedrals_bis(2,di), os(oi).dihedrals_ibs(2,di) ) );
            if iu1 == 0 || iu2 == 0
                error("Dihedral bead does not exist.")
            end
            bond(1) = bond_type;
            bond(2) = iu1;
            bond(3) = iu2;
            if ~any(all(bonds == bond, 1))
                bond_count = bond_count + 1;
                bonds(:,bond_count) = bond;
            end
            iu3 = get_iu( oi, os(oi).get_io( os(oi).dihedrals_bis(3,di), os(oi).dihedrals_ibs(3,di) ) );
            iu4 = get_iu( oi, os(oi).get_io( os(oi).dihedrals_bis(4,di), os(oi).dihedrals_ibs(4,di) ) );
            if iu3 == 0 || iu4 == 0
                error("Dihedral bead does not exist.")
            end
            bond(1) = bond_type;
            bond(2) = iu3;
            bond(3) = iu4;
            if ~any(all(bonds == bond, 1))
                bond_count = bond_count + 1;
                bonds(:,bond_count) = bond;
            end
        end
    end

    %%% rigid bonds for reactions
    for ri = 1:length(rs)
        nbond = size(rs(ri).bonds_init,1);
        for bondi = 1:nbond
            bond_type = rs(ri).bonds_init(bondi,1);
            ir1 = rs(ri).bonds_init(bondi,2);
            ir2 = rs(ri).bonds_init(bondi,3);
            if ir1 <= rs(ri).n_site5
                for si = 1:rs(ri).nsite5
                    oi = rs(ri).site5s_oi(si);
                    bi = rs(ri).site5s_bi(si);
                    ib = rs(ri).site5s_ibs(ir1,si);
                    io = os(oi).get_io(bi,ib);
                    iu1 = get_iu(oi,io);
                    oi = rs(ri).site5s_oi(si);
                    bi = rs(ri).site5s_bi(si);
                    ib = rs(ri).site5s_ibs(ir2,si);
                    io = os(oi).get_io(bi,ib);
                    iu2 = get_iu(oi,io);
                    bond_count = bond_count + 1;
                    bonds(1,bond_count) = bond_type;
                    bonds(2,bond_count) = iu1;
                    bonds(3,bond_count) = iu2;
                end
            else
                ir1 = ir1 - rs(ri).n_site5;
                ir2 = ir2 - rs(ri).n_site5;
                for si = 1:rs(ri).nsite3
                    oi = rs(ri).site3s_oi(si);
                    bi = rs(ri).site3s_bi(si);
                    ib = rs(ri).site3s_ibs(ir1,si);
                    io = os(oi).get_io(bi,ib);
                    iu1 = get_iu(oi,io);
                    oi = rs(ri).site3s_oi(si);
                    bi = rs(ri).site3s_bi(si);
                    ib = rs(ri).site3s_ibs(ir2,si);
                    io = os(oi).get_io(bi,ib);
                    iu2 = get_iu(oi,io);
                    bond_count = bond_count + 1;
                    bonds(1,bond_count) = bond_type;
                    bonds(2,bond_count) = iu1;
                    bonds(3,bond_count) = iu2;
                end
            end
        end
    end

    %%% compile angle info
    angles = zeros(4,0);
    angle_count = 0;
    for oi = 1:length(os)
        for ai = 1:length(os(oi).angles_apot)
            if os(oi).angles_status(ai) == 1
                angle_type = os(oi).angles_apot(ai).index;
                iu1 = get_iu( oi, os(oi).get_io( os(oi).angles_bis(1,ai), os(oi).angles_ibs(1,ai) ) );
                iu2 = get_iu( oi, os(oi).get_io( os(oi).angles_bis(2,ai), os(oi).angles_ibs(2,ai) ) );
                iu3 = get_iu( oi, os(oi).get_io( os(oi).angles_bis(3,ai), os(oi).angles_ibs(3,ai) ) );
                iu4 = get_iu( oi, os(oi).get_io( os(oi).angles_bis(4,ai), os(oi).angles_ibs(4,ai) ) );
                angle_count = angle_count + 1;
                angles(1,angle_count) = angle_type;
                angles(2,angle_count) = iu1;
                angles(3,angle_count) = iu2;
                angles(4,angle_count) = iu3;
                angle_count = angle_count + 1;
                angles(1,angle_count) = angle_type;
                angles(2,angle_count) = iu2;
                angles(3,angle_count) = iu3;
                angles(4,angle_count) = iu4;
            end
        end
    end

    %%% compile dihedral info
    dihedrals = zeros(5,0);
    dihedral_count = 0;
    for oi = 1:length(os)
        for di = 1:length(os(oi).dihedrals_dpot)
            dihedral_type = os(oi).dihedrals_dpot(di).index;
            iu1 = get_iu( oi, os(oi).get_io( os(oi).dihedrals_bis(1,di), os(oi).dihedrals_ibs(1,di) ) );
            iu2 = get_iu( oi, os(oi).get_io( os(oi).dihedrals_bis(2,di), os(oi).dihedrals_ibs(2,di) ) );
            iu3 = get_iu( oi, os(oi).get_io( os(oi).dihedrals_bis(3,di), os(oi).dihedrals_ibs(3,di) ) );
            iu4 = get_iu( oi, os(oi).get_io( os(oi).dihedrals_bis(4,di), os(oi).dihedrals_ibs(4,di) ) );
            dihedral_count = dihedral_count + 1;
            dihedrals(1,dihedral_count) = dihedral_type;
            dihedrals(2,dihedral_count) = iu1;
            dihedrals(3,dihedral_count) = iu2;
            dihedrals(4,dihedral_count) = iu3;
            dihedrals(5,dihedral_count) = iu4;
        end
    end

    %%% write simulation geometry file
    charges = zeros(1,atom_count);
    ars.writeGeo(geoFile,dbox,atoms,bonds,angles,dihedrals=dihedrals,masses=masses,charges=charges,nbondType=nABADtype(2),nangleType=nABADtype(3),ndihedralType=nABADtype(4));

    %%% write visualization geometry file
    ars.writeGeo(geoVisFile,dbox,atoms_vis,bonds,angles,dihedrals=dihedrals)
end


%%% write lammps input file
function write_input(inputFile,reactFold,p,nABADtype,ls,rs,pots,apots,dpots)
    disp("Writing input file...")

    %%% calculate communication cutoff
    if p.comm_cutoff == 0
        U_max = 12;
        p.comm_cutoff = p.r12_cut_WCA;
        if numEntries(pots) > 0
            for p_name = keys(pots)'
                r12_max = pots(p_name).calc_separation(U_max);
                p.comm_cutoff = max([p.comm_cutoff,r12_max]);
            end
        end
        for li = 1:length(ls)
            r12_cut = ls(li).r12_cut;
            p.comm_cutoff = max([p.comm_cutoff,r12_cut]);
        end
        for ri = 1:length(rs)
            r12_min_max = max(rs(ri).r12s_min);
            r12_max_max = max(rs(ri).r12s_max);
            p.comm_cutoff = max([p.comm_cutoff,r12_min_max,r12_max_max]);
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
        "atom_style      full\n",...
        "read_data       geometry.in &\n",...
        "                extra/bond/per/atom 10 &\n",...
        "                extra/angle/per/atom 10 &\n",...
        "                extra/dihedral/per/atom 10 &\n",...
        "                extra/special/per/atom 10\n\n"));

    %%% neighbor list
    fprintf(f,strcat(...
        "## Parameters\n",...
        "neighbor        ", ars.fstring(p.verlet_skin,0,2), " bin\n",...
        "neigh_modify    every ", num2str(p.neigh_every), "\n",...
        "neigh_modify    exclude molecule/intra all\n"));

    %%% pairwise interactions
    fprintf(f,strcat(...
        "pair_style      hybrid/overlay lj/cut ", ars.fstring(p.r12_cut_WCA,0,2), " zero ", ars.fstring(p.comm_cutoff,0,2), "\n",...
        "pair_coeff      * * zero\n",...
        "pair_coeff      1 1 lj/cut ", ars.fstring(p.epsilon,0,2), " ", ars.fstring(p.sigma,0,2), "\n",...
        "pair_modify     pair lj/cut shift yes\n"));

    %%% bonds
    if nABADtype(2) > 0
        bond_styles = unique([values(pots).style]);
        if isscalar(bond_styles)
            is_hybrid = false;
            fprintf(f,strcat(...
                "bond_style      ", bond_styles{1}, "\n"));
        else
            is_hybrid = true;
            fprintf(f,...
                "bond_style      hybrid");
            for si = 1:length(bond_styles)
                fprintf(f,strcat(" ", bond_styles{si}));
            end
            fprintf(f,"\n");
        end
        for p_name = keys(pots)'
            pots(p_name).write_potential(f,is_hybrid);
        end
    end

    %%% angles
    if nABADtype(3) > 0
        angle_styles = unique([values(apots).style]);
        if isscalar(angle_styles)
            is_hybrid = false;
            fprintf(f,strcat(...
                "angle_style     ", angle_styles{1}, "\n"));
        else
            is_hybrid = true;
            fprintf(f,...
                "angle_style     hybrid");
            for si = 1:length(angle_styles)
                fprintf(f,strcat(" ", angle_styles{si}));
            end
            fprintf(f,"\n");
        end
        for a_name = keys(apots)'
            apots(a_name).write_potential(f,is_hybrid);
        end
    end

    %%% dihedrals
    if nABADtype(4) > 0
        dihedral_styles = unique([values(dpots).style]);
        if isscalar(dihedral_styles)
            is_hybrid = false;
            fprintf(f,strcat(...
                "dihedral_style  ", dihedral_styles{1}, "\n"));
        else
            is_hybrid = true;
            fprintf(f,...
                "dihedral_style  hybrid");
            for si = 1:length(dihedral_styles)
                fprintf(f,strcat(" ", dihedral_styles{si}));
            end
            fprintf(f,"\n");
        end
        for d_name = keys(dpots)'
            dpots(d_name).write_potential(f,is_hybrid);
        end
    end

    %%% charge varaible
    fprintf(f,...
        "variable        varQ atom q\n");

    %%% patch group
    if nABADtype(1) >= 2
        fprintf(f,strcat(...
            "group           patchy type 2:", num2str(nABADtype(1)), "\n"));
    end

    %%% sticky patch group
    if nABADtype(1) >= 3
        fprintf(f,strcat(...
            "group           sticky type 3:", num2str(nABADtype(1)), "\n"));
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

    %%% linker reactions
    if ~isempty(ls)
        fprintf(f,"## Reactions\n");
        for li = 1:length(ls)
            ls(li).write_bond_create(f,p.react_every);
        end
        fprintf(f,"\n");
    end

    %%% reacter reactions
    if ~isempty(rs)
        fprintf(f,"## Reactions\n");
        ars.createEmptyFold(reactFold);
        for ri = 1:length(rs)
            rs(ri).write_react_files(reactFold,pots,apots,dpots)
            rs(ri).write_molecules(f)
        end
        fprintf(f,"fix             reactions sticky bond/react reset_mol_ids no");
        for ri = 1:length(rs)
            rs(ri).write_react(f,p.comm_cutoff,p.react_every)
        end
        fprintf(f,"\n\n");
    end
    
    %%% updates
    if p.nstep > 0
        fprintf(f,strcat(...
            "## Updates\n",...
            "dump            dumpT all custom ", num2str(p.dump_every), " trajectory.dat id mol xs ys zs\n",...
            "dump_modify     dumpT sort id\n",...
            "compute         compB1 sticky bond/local dist engpot\n",...
            "compute         compB2 sticky property/local btype batom1 batom2\n",...
            "dump            dumpB sticky local ", num2str(p.dump_every), " dump_bonds.dat index c_compB1[1] c_compB1[2] c_compB2[1] c_compB2[2] c_compB2[3]\n"));
        if nABADtype(3) > 0
            fprintf(f,strcat(...
			"compute         compA1 patchy angle/local theta eng\n",...
			"compute         compA2 patchy property/local atype aatom1 aatom2 aatom3\n",...
		    "dump            dumpA patchy local ", num2str(p.dump_every), " dump_angles.dat index c_compA1[1] c_compA1[2] c_compA2[1] c_compA2[2] c_compA2[3] c_compA2[4]\n"));
        end
        if nABADtype(4) > 0
            fprintf(f,strcat(...
			"compute         compD1 patchy dihedral/local phi\n",...
			"compute         compD2 patchy property/local dtype datom1 datom2 datom3 datom4\n",...
		    "dump            dumpD patchy local ", num2str(p.dump_every), " dump_dihedrals.dat index c_compD1 c_compD2[1] c_compD2[2] c_compD2[3] c_compD2[4] c_compD2[5]\n"));
        end
        if any([rs.is_charged])
        fprintf(f,strcat(...
            "dump            dumpC sticky custom ", num2str(p.dump_every), " dump_charges.dat id q\n",...
            "dump_modify     dumpC sort id\n"));
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


%%% calculate origami indices of origami type
function ois = get_ois(o_counts,o_name)
    o_count = 0;
    for o_name_key = keys(o_counts)'
        if o_name_key == o_name
            ois = o_count+1:o_count+o_counts(o_name);
            return
        end
        o_count = o_count + o_counts(o_name_key);
    end
end