%%% Housekeeping
clc; clear; close all;
rng(43)

%%% To Do
% replace bond write/break with react if unlinking is desired.
% add option to set linker and connection spring constant.
% add option to initialize block from two flexible connections.

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
% connection - permenant (usually ssDNA scaffold) bond
% linker - switchable (usually hybridizing DNA) bond

%%% Input file
% blocks: rigid bodies
  % name - how to identify the block for patches and origamis
  % pattern - arangement of beads in the xy-plane (see block object)
  % height - number of beads in the z-direction
% patches: locations on blocks used for bonded interactions
  % name - how to identify the patch for connections and linkers
  % block - name of block on which to place the patch
  % theta - polar angle in xy-plane (in degrees) of patch location
  % radius - distance from z-axis (in r12_helix units) of patch location
  % z - height along z-axis (in r12_bead units) of patch location
% origami: collection of connected blocks
  % name - how to identify the origami for defining parameters and adding linkers
  % block - names of block types (blocks are indexed in the given order)
  % conn - permanent bond between indexed blocks at given location
  % count - number of origamis to create
% linker: breakable bond
  % origami - name of origamis on which to create the linker
  % block indices - index, or "A" for all blocks, or "B" for all but the last block
  % location - where to place the ends of the linker
% note: locations are defined by a flag (B or P, for bead or patch),
  % followed by a hyphen, followed by a bead or patch identifier; for
  % beads, the identifier is the helix (L/M/R for left/middle/right)
  % directly followed by the height; for patches the identifier is the name
  % of the patch.
% note: linkers require two lines, one starting with the keyword "linker"
  % that defines the location of the 5' end, and another line directly
  % afterwards that defines the location of the 3' end.


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% read input
inFile = "./designs/control.txt";
[os,origami_types,linker_types,potentials] = read_input(inFile);

%%% output parameters
outFold = "/Users/dduke/Files/block_tether/free_origami/lammps/active/";
nsim = 1;

%%% simulation parameters
nstep_eq        = 1E4;      % steps         - if and how long to equilibrate/shrink
shrink_ratio    = 1;        % none          - box compression (final/initial)
nstep_prod      = 3E7;      % steps         - if and how long to run producton
dump_every      = 1E5;      % steps         - how often to write to output

%%% computational parameters
dt              = 0.08;     % ns            - time step
dbox            = 320;      % nm            - periodic boundary diameter
verlet_skin     = 4;        % nm            - width of neighbor list skin
neigh_every     = 1E1;      % steps         - how often to consider updating neighbor list
react_every     = 1E1;      % steps         - how often to check for linker hybridization

%%% physical parameters
T               = 300;      % K             - temperature
sigma           = 10;       % nm            - WCA distance parameter
epsilon         = 1;        % kcal/mol      - WCA energy parameter
r12_bead        = 5;        % nm            - helix separation
r12_helix       = 5;        % nm            - bead separation
k_x_conn        = 1;        % kcal/mol/nm2  - connection spring constant
k_x_linker      = 1;        % kcal/mol/nm2  - linker spring constant
U_cut_linker    = 8;        % kcal/mol      - linker reaction cutoff
k_theta         = 6;        % kcal/mol/rad2 - angle spring constant

%%% create parameters class
p = parameters(nstep_eq,shrink_ratio,nstep_prod,dump_every,...
               dt,dbox,verlet_skin,neigh_every,react_every,...
               T,sigma,epsilon,r12_bead,r12_helix,...
               k_x_conn,k_x_linker,U_cut_linker,k_theta);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Write LAMMPS Files %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    compose_geo(geoFile,geoVisFile,os,origami_types,linker_types,dbox);
    
    %%% write lammps input file
    inputFile = simFold + "lammps.in";
    write_input(inputFile,p,origami_types,linker_types,potentials)
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Simulation Functions %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% initilize positions of the entire system
function os = init_positions(os,p)
    max_attempts = 10;

    %%% initialize avoided positions
    r_other = [];
    r_other_initial = r_other;

    %%% system attempt loop
    attempts = 0;
    while true

        %%% check attempts
        if attempts == max_attempts
            error("Could not place origamis.")
        end

        %%% loop over origamis
        for oi = 1:length(os)

            %%% add origami
            disp(strcat("Initializing origami ",num2str(oi),"..."))
            [os(oi),failed,r_other] = os(oi).init_positions(p,r_other);

            %%% reset if failed
            if failed == true
                attempts = attempts + 1;
                r_other = r_other_initial;
                break
            end
        end

        %%% reset if failed
        if failed == true
            continue
        end

        %%% system successfully initiated
        break
    end

    %%% report the wonderful news
    fprintf("Initialization complete.\n")
end


%%% map between universal index and origami indices
function [get_oi,get_io,get_iu] = map_indices(os)
    max_io = max([os.n]);
    n_uni = sum([os.n]);
    get_oi = zeros(1,n_uni);
    get_io = zeros(1,n_uni);
    get_iu = zeros(length(os),max_io);
    iu = 0;
    for oi = 1:length(os)
        for io = 1:os(oi).n
            iu = iu+1;
            get_oi(iu) = oi;
            get_io(iu) = io;
            get_iu(oi,io) = iu;
        end
    end
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% File Handlers %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% read input file and create corresponding origami objects
function [os,origami_types,linker_types,potential_types] = read_input(inFile)

    %%% open file
    f = fopen(inFile, 'r');
    if f == -1
        error("Could not open file.");
    end

    %%% initialize dictionaries
    potential_types = dictionary();
    block_templates = dictionary();
    origami_types = dictionary();
    origami_templates = dictionary();
    linker_types = dictionary();

    %%% loop over lines
    while ~feof(f)
        line = fgetl(f);

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

        %%% set the appropriate value
        switch extract{1}
            case 'potential'
                potential.type = extract{3};
                potential.r12_eq = str2double(extract{4});
                potential.parameters = str2double(extract(5:end));
                potential_types(extract{2}) = potential;
            case 'block'
                block_templates(extract{2}) = block(convertCharsToStrings(extract{3}),str2double(extract{4}));
            case 'patch'
                block_templates(extract{3}) = block_templates(extract{3}).add_patch(extract{2},str2double(extract{4}),str2double(extract{5}),str2double(extract{6}));
            case 'origami'
                origami_templates(extract{2}) = origami(string(extract{2}));
                origami_type.count = str2double(extract{3});
                origami_type.conns_r12_eq = [];
                origami_type.angles_theta_eq = [];
                origami_types(extract{2}) = origami_type;
            case 'linker'
                linker_type.r12_eq = str2double(extract{3});
                linker_types(extract{2}) = linker_type;
            otherwise
                if isKey(origami_types,extract{1})
                    switch extract{2}
                        case 'blocks'
                            for bi = 1:length(extract)-2
                                b = block_templates(extract{bi+2});
                                origami_templates(extract{1}) = origami_templates(extract{1}).add_block(b);
                            end
                        case 'rigid'
                            origami_templates(extract{1}).rigid = string(extract{3});
                        case 'conn'
                            origami_templates(extract{1}) = origami_templates(extract{1}).add_conn(str2double(extract{3}),extract{4},str2double(extract{5}),extract{6},str2double(extract{7}));
                            origami_types(extract{1}).conns_r12_eq = [ origami_types(extract{1}).conns_r12_eq potential_types(extract{7}).r12_eq ];
                        case 'angle'
                            theta_init = str2double(extract{11});
                            if length(extract) > 11
                                theta_init = str2double(extract{12});
                            end
                            origami_templates(extract{1}) = origami_templates(extract{1}).add_angle(str2double(extract{3}),extract{4},str2double(extract{5}),extract{6},str2double(extract{7}),extract{8},str2double(extract{9}),extract{10},str2double(extract{11}),theta_init);
                            origami_types(extract{1}).angles_theta_eq = [ origami_types(extract{1}).angles_theta_eq str2double(extract{11}) ];
                        otherwise
                            error("Unknown origami parameter: " + extract{2})
                    end
                elseif isKey(linker_types,extract{1})
                    switch extract{2}
                        case '5p'
                            origami_templates(extract{3}) = origami_templates(extract{3}).add_linker(string(extract{1}),1,extract{4},extract{5});
                        case '3p'
                            origami_templates(extract{3}) = origami_templates(extract{3}).add_linker(string(extract{1}),0,extract{4},extract{5});
                        otherwise
                            error("Unknown linker parameter: " + extract{2})
                    end
                else
                    error("Unknown system parameter: " + extract{1})
                end
        end
    end
    fclose(f);

    %%% create origamis
    os = origami.empty;
    for o_name = keys(origami_types)'
        count = origami_types(o_name).count;
        os(length(os)+1:length(os)+count) = origami_templates(o_name);
    end
end


%%% write lammps geometry file
function compose_geo(geoFile,geoVisFile,os,origami_types,linker_types,dbox)

    %%% grouping atoms
    % atom ID - universal bead index
    % molecule tag - universal block index
    % atom type - unlinked bead (1), unlinked patch (2), linker bead/patch (3:2+nlinker*2)
    % bond type - linker (1:nlinker), connection (nlinker+1:nlinker+nconn)

    %%% calculate index mapping
    [get_oi,get_io,get_iu] = map_indices(os);

    %%% get counts
    natom = sum([os.n]);
    nconn = sum([os.nconn]);
    nangle = sum([os.nangle]);

    %%% count bond types
    nlinker = numEntries(linker_types);
    nconnType = length([values(origami_types).conns_r12_eq]);
    nbondType = nlinker+nconnType+1;

    %%% initialize mass info
    natomType = 2 + nlinker*2;
    masses = ones(1,natomType);
    mass_patch = 0.01;
    masses(2) = mass_patch;

    %%% compile atom info
    atoms = zeros(5,natom);
    rigid_count = 0;
    for iu = 1:natom
        oi = get_oi(iu);
        io = get_io(iu);
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
        atoms(1,iu) = rigid_count;
        if os(oi).is_patch(io)
            atoms(2,iu) = 2;
        else
            atoms(2,iu) = 1;
        end
        atoms(3:5,iu) = os(oi).r(:,io);
    end
    for oi = 1:length(os)
        for lio = 1:os(oi).nlink5
            io = os(oi).link5s_io(lio);
            iu = get_iu(oi,io);
            o_name = os(oi).link5s_name(lio);
            li = find(keys(linker_types)==o_name);
            atomType = 1 + 2*li;
            atoms(2,iu) = atomType;
            if os(oi).is_patch(io)
                masses(atomType) = mass_patch;
            end
        end
        for lio = 1:os(oi).nlink3
            io = os(oi).link3s_io(lio);
            iu = get_iu(oi,io);
            o_name = os(oi).link3s_name(lio);
            li = find(keys(linker_types)==o_name);
            atomType = 2 + 2*li;
            atoms(2,iu) = atomType;
            if os(oi).is_patch(io)
                masses(atomType) = mass_patch;
            end
        end
    end

    %%% connection bonds
    bonds = zeros(3,0);
    if nconn > 0
        bond_type = nlinker;
        bond_count = 0;
        for o_name = keys(origami_types)'
            ois = find([os.name]==o_name);
            for ci = 1:length(origami_types(o_name).conns_r12_eq)
                bond_type = bond_type + 1;
                for oi = ois
                    bond_count = bond_count + 1;
                    iu1 = get_iu( oi, os(oi).get_io( os(oi).conns_bis(1,ci), os(oi).conns_ibs(1,ci) ) );
                    iu2 = get_iu( oi, os(oi).get_io( os(oi).conns_bis(2,ci), os(oi).conns_ibs(2,ci) ) );
                    bonds(1,bond_count) = bond_type;
                    bonds(2,bond_count) = iu1;
                    bonds(3,bond_count) = iu2;
                end
            end
        end
    end

    %%% rigid bonds for angles
    if nangle > 0
        bond_type = nbondType;
        for o_name = keys(origami_types)'
            ois = find([os.name]==o_name);
            for ai = 1:length(origami_types(o_name).angles_theta_eq)
                for oi = ois
                    iu1 = get_iu( oi, os(oi).get_io( os(oi).angles_bis(1,ai), os(oi).angles_ibs(1,ai) ) );
                    iu2 = get_iu( oi, os(oi).get_io( os(oi).angles_bis(2,ai), os(oi).angles_ibs(2,ai) ) );
                    bond_count = bond_count + 1;
                    bonds(1,bond_count) = bond_type;
                    bonds(2,bond_count) = iu1;
                    bonds(3,bond_count) = iu2;
                    iu3 = get_iu( oi, os(oi).get_io( os(oi).angles_bis(3,ai), os(oi).angles_ibs(3,ai) ) );
                    iu4 = get_iu( oi, os(oi).get_io( os(oi).angles_bis(4,ai), os(oi).angles_ibs(4,ai) ) );
                    bond_count = bond_count + 1;
                    bonds(1,bond_count) = bond_type;
                    bonds(2,bond_count) = iu3;
                    bonds(3,bond_count) = iu4;
                end
            end
        end
    end

    %%% compile angle info
    angles = zeros(4,nangle);
    angle_type = 0;
    if nangle > 0
        angle_count = 0;
        for o_name = keys(origami_types)'
            ois = find([os.name]==o_name);
            for ai = 1:length(origami_types(o_name).angles_theta_eq)
                angle_type = angle_type + 1;
                for oi = ois
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
    end
    nangleType = angle_type+1;

    %%% write simulation geometry file
    ars.writeGeo(geoFile,dbox,atoms,bonds,angles,masses=masses,nbondType=nbondType,nangleType=nangleType)

    %%% compile atom info for visualization
    atoms_vis = atoms;
    norigami_type = numEntries(origami_types);
    for iu = 1:natom
        oi = get_oi(iu);
        io = get_io(iu);
        atoms_vis(1,iu) = oi;
        if os(oi).is_patch(io)
            atoms_vis(2,iu) = norigami_type + 1;
        else
            o_name = os(oi).name;
            atoms_vis(2,iu) = find(keys(origami_types)==o_name);
        end
    end
    for oi = 1:length(os)
        for lio = 1:os(oi).nlink5
            io = os(oi).link5s_io(lio);
            iu = get_iu(oi,io);
            li = find(keys(linker_types)==os(oi).link5s_name(lio));
            atoms_vis(2,iu) = norigami_type + 1 + li;
        end
        for lio = 1:os(oi).nlink3
            io = os(oi).link3s_io(lio);
            iu = get_iu(oi,io);
            li = find(keys(linker_types)==os(oi).link3s_name(lio));
            atoms_vis(2,iu) = norigami_type + 1 + li;
        end
    end

    %%% write visualization geometry file
    ars.writeGeo(geoVisFile,dbox,atoms_vis,bonds,angles)
end


%%% write lammps input file
function write_input(inputFile,p,origami_types,linker_types,potential_types)

    %%% get counts
    nlinker = numEntries(linker_types);
    natomType = 2+nlinker*2;
    nconnType = length([values(origami_types).conns_r12_eq]);
    nbondType = nlinker+nconnType+1;
    nangleType = length([values(origami_types).angles_theta_eq]);

    %%% calculate communication cutoff
    U_max = 12;
    comm_cutoff = p.r12_cut_WCA;
    for o_name = keys(origami_types)'
        for r12_eq = origami_types(o_name).conns_r12_eq
            r12_max = r12_eq + sqrt(2*U_max/p.k_x_conn);
            comm_cutoff = max([comm_cutoff,r12_max]);
        end
    end
    if nlinker > 0
        for l_name = keys(linker_types)'
            r12_cut = linker_types(l_name).r12_eq + sqrt(2*p.U_cut_linker/p.k_x_linker);
            r12_max = linker_types(l_name).r12_eq + sqrt(2*U_max/p.k_x_linker);
            comm_cutoff = max([comm_cutoff,r12_cut,r12_max]);
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

    %%% prepare bonds
    if nbondType > 1
        fprintf(f,...
            "bond_style      hybrid harmonic zero\n");
        bond_type = 0;
    else
        fprintf(f,...
            "bond_style      zero\n");
    end


    %%% linker bonds
    if nlinker > 0
        for l_name = keys(linker_types)'
            bond_type = bond_type + 1;
            r12_eq = linker_types(l_name).r12_eq;
            fprintf(f,strcat(...
                "bond_coeff      ", num2str(bond_type), " harmonic ", ars.fstring(p.k_x_linker/2,0,2), " ", ars.fstring(r12_eq,0,2), "\n"));
        end
    end

    %%% connection bonds
    if nconnType > 0
        for o_name = keys(origami_types)'
            for r12_eq = origami_types(o_name).conns_r12_eq
                bond_type = bond_type + 1;
                fprintf(f,strcat(...
                    "bond_coeff      ", num2str(bond_type), " harmonic ", ars.fstring(p.k_x_conn/2,0,2), " ", ars.fstring(r12_eq,0,2), "\n"));
            end
        end
    end

    %%% dummy bond
    if nbondType > 1
        bond_type = bond_type + 1;
        fprintf(f,strcat(...
            "bond_coeff      ", num2str(bond_type), " zero\n"));
    else
        fprintf(f,strcat(...
            "bond_coeff      1\n"));
    end

    %%% angles
    if nangleType > 0
        fprintf(f,...
            "angle_style     harmonic\n");
        angle_type = 0;
        for o_name = keys(origami_types)'
            for theta_eq = origami_types(o_name).angles_theta_eq
                angle_type = angle_type + 1;
                fprintf(f,strcat(...
                    "angle_coeff     ", num2str(angle_type), " ", ars.fstring(p.k_theta/2,0,2), " ", ars.fstring(theta_eq,0,2), "\n"));
            end
        end
    end

    %%% linker angle
    if nlinker > 0
        if nangleType == 0
            angle_type = 1;
        else
            angle_type = angle_type + 1;
        end
        fprintf(f,strcat(...
            "angle_coeff     ", num2str(angle_type), " ", ars.fstring(100/2,0,2), " ", ars.fstring(180,0,2), "\n"));
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
    if p.nstep_eq > 0
        fprintf(f,strcat(...
            "## Relaxation\n",...
            "timestep        ", num2str(p.dt/10), "\n"));
        if p.shrink_ratio ~= 1
            fprintf(f,strcat(...
                "fix             deformation all deform 1 &\n",...
                "                x final ", ars.fstring(-p.dbox/2*p.shrink_ratio,0,2), " ", ars.fstring(p.dbox/2*p.shrink_ratio,0,2), " &\n",...
                "                y final ", ars.fstring(-p.dbox/2*p.shrink_ratio,0,2), " ", ars.fstring(p.dbox/2*p.shrink_ratio,0,2), " &\n",...
                "                z final ", ars.fstring(-p.dbox/2*p.shrink_ratio,0,2), " ", ars.fstring(p.dbox/2*p.shrink_ratio,0,2), "\n"));
        end
        fprintf(f,strcat(...
            "run             ", num2str(p.nstep_eq), "\n"));
        if p.shrink_ratio ~= 1
            fprintf(f,strcat(...
                "unfix           deformation\n"));
        end
        fprintf(f,strcat(...
            "reset_timestep  0\n\n"));
    end

    %%% reactions
    if nlinker > 0
        fprintf(f,"## Reactions\n");
        for l_name = keys(linker_types)'
            li = find(keys(linker_types)==l_name);
            r12_cut = linker_types(l_name).r12_eq + sqrt(2*p.U_cut_linker/p.k_x_linker);
            fixName = strcat("linker",num2str(li));
            add_bond_create(f,fixName,"linker",p.react_every,li*2+1,li*2+2,r12_cut,li,angle_type);
        end
        fprintf(f,"\n");
    end
    
    %%% updates
    if p.nstep_prod > 0
        fprintf(f,strcat(...
            "## Updates\n",...
            "dump            dump1 all custom ", num2str(p.dump_every), " trajectory.dat id mol xs ys zs\n",...
            "dump_modify     dump1 sort id\n",...
            "compute         compD1a all bond/local dist engpot\n",...
            "compute         compD1b all property/local btype batom1 batom2\n",...
            "dump            dumpD1 all local ", num2str(p.dump_every), " dump_bonds.txt index c_compD1a[1] c_compD1a[2] c_compD1b[1] c_compD1b[2] c_compD1b[3]\n",...
			"compute         compD2a all angle/local theta eng\n",...
			"compute         compD2b all property/local atype aatom1 aatom2 aatom3\n",...
		    "dump            dumpD2 all local ", num2str(p.dump_every), " dump_angles.dat index c_compD2a[1] c_compD2a[2] c_compD2b[1] c_compD2b[2] c_compD2b[3] c_compD2b[4]\n\n"));
    end

    %%% production
    if p.nstep_prod > 0
        fprintf(f,strcat(...
            "## Production\n",...
            "timestep        ", num2str(p.dt), "\n",...
            "run             ", num2str(p.nstep_prod), "\n"));
    end
    
    %%% finalize and close file
    fprintf(f,"\n#------- End Input -------#\n\n");
    fclose(f);
end


%%% write fix bond create command
function add_bond_create(f,fixName,groupName,react_every,atomType1,atomType2,r12_cut,bondType,angleType)
    fprintf(f,strcat(...
        "fix             ", fixName, " ", groupName, " ",...
        "bond/create/angle ", num2str(react_every), " ",...
        num2str(atomType1), " ",num2str(atomType2), " ",...
        ars.fstring(r12_cut,0,2), " ", num2str(bondType), " "));
    if angleType ~= 0 
        fprintf(f,strcat("atype ", num2str(angleType), " aconstrain 120 180 "));
    end
    fprintf(f,strcat(...
        "iparam 1 ", num2str(atomType1), " ",...
        "jparam 1 ", num2str(atomType2), "\n"));
end
