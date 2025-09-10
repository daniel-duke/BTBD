%%% Housekeeping
clc; clear; close all;
rng(42)

%%% To Do
% correct pairwise coffecients for linked block beads.
% replace bond write/break with react if unlinking is desired.
% add option to define origami initial configuration.
% add option to set linker spring constnat.
% add option to set connection distance and spring constant.
% initialize origami in empty box (ensuring no internal overlap), then
  % place into system (ensuring no overlap with other origamis).
% improve bead notation: bi = block index, ib = bead index within block,
  % thus removing need for confusing pi and ui indices.

%%% Notation
% r - position vector
% r12 - vector pointing from position of bead 1 to bead 2
% o, os, oi - origami, list of origamis, origami index
% b, bs, bi - block, list of blocks, block index
% nvar - number of var
% var.n or n_var - number of beads in var
% i - bead index within block
% pi - bead (particle) index across entire origami
% ui - (universal) bead index across entire system
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
  % radius - distance from z-axis (in r12_adj_block units) of patch location
  % z - height along z-axis (in r12_eq_block units) of patch location
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
inFile = "./designs/triarm.txt";
[os,linked,o_types,conn_types,conn_r12_eqs] = read_input(inFile);

%%% output parameters
outFold = "/Users/dduke/Files/block_tether/network/active/";
nsim = 1;

%%% simulation parameters
nstep_eq        = 1E4;      % steps         - if and how long to equilibrate/shrink
shrink_ratio    = 1;        % none          - box compression (final/initial)
nstep_prod      = 1E7;      % steps         - if and how long to run producton
dump_every      = 1E4;      % steps         - how often to write to output

%%% computational parameters
dt              = 0.01;     % ns            - time step
dbox            = 150;      % nm            - periodic boundary diameter
verlet_skin     = 4;        % nm            - width of neighbor list skin (= r12_cut - r12_cut_WCA)
neigh_every     = 1E1;      % steps         - how often to consider updating neighbor list
react_every     = 1E1;      % steps         - how often to check for linker hybridization

%%% physical parameters
T               = 300;      % K             - temperature
sigma           = 10;       % nm            - WCA distance parameter
epsilon         = 0.1;      % kcal/mol      - WCA energy parameter
r12_eq_block    = 5;        % nm            - equilibrium block bead separation
r12_adj_block   = 5;        % nm            - adjacent block helix separation
k_x_conn        = 0.1;        % kcal/mol/nm2  - connection spring constant
r12_eq_linker   = 3;        % nm            - linker equilibrium distance
k_x_linker      = 1;        % kcal/mol/nm2  - linker spring constant

%%% create parameters class
p = parameters(nstep_eq,shrink_ratio,nstep_prod,dump_every,...
               dt,dbox,verlet_skin,neigh_every,react_every,...
               T,sigma,epsilon,r12_eq_block,r12_adj_block,...
               k_x_conn,r12_eq_linker,k_x_linker);


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
    os = place_origami(os,p);

    %%% write lammps simulation geometry file
    geoFile = simFold + "geometry.in";
    geoVisFile = simFold + "geometry_vis.in";
    compose_geo(geoFile,geoVisFile,os,linked,o_types,conn_types,dbox);
    
    %%% write lammps input file
    inputFile = simFold + "lammps.in";
    write_input(inputFile,p,linked,conn_types,conn_r12_eqs)
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Simulation Functions %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% initilize positions of the entire system
function os = place_origami(os,p)
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
            [os(oi),fail,r_other] = os(oi).init(p,r_other);

            %%% reset if failed
            if fail == true
                attempts = attempts + 1;
                r_other = r_other_initial;
                break
            end

            %%% update positions
            os(oi).r = os(oi).update_r;
        end

        %%% reset if failed
        if fail == true
            continue
        end

        %%% system successfully initiated
        break
    end

    %%% report the wonderful news
    fprintf("Initialization complete.\n")
end


%%% create linkers by updating linked matrix
function linked = linker(oi1s,bi1s,loc1,oi2s,bi2s,loc2,os,linked)

    %%% meaning of linked matrix values
    % 0 - not linked to anything
    % odd - linked to li+1
    % even - linked to li-1

    %%% identify if linker(s) are intra or inter
    if isequal(oi1s,oi2s)
        interlink = false;
        linker_index = ars.myMax(linked)/2;
    else
        interlink = true;
        linker_index = ars.myMax(linked)/2 + 1;
    end

    %%% set values in linker array
    for oi1 = oi1s
        for oi2 = oi2s
            if oi2 == oi1 || interlink == true
                if interlink == false
                    linker_index = linker_index + 1;
                end
                for j1 = 1:length(bi1s)
                    for j2 = 1:length(bi2s)
                        bi1 = bi1s(j1);
                        bi2 = bi2s(j2);
                        i1 = os(oi1).bs(bi1).interpret_loc(loc1);
                        i2 = os(oi2).bs(bi2).interpret_loc(loc2);
                        pi1 = os(oi1).get_pi(bi1,i1);
                        pi2 = os(oi2).get_pi(bi2,i2);
                        linked(oi1,pi1) = linker_index*2-1;
                        linked(oi2,pi2) = linker_index*2;
                    end
                end
            end
        end
    end
end


%%% map between universal index and origami/particle index
function [get_oi,get_pi,get_ui] = map_indices(os)
    max_pi = max([os.n]);
    n_uni = sum([os.n]);
    get_oi = zeros(1,n_uni);
    get_pi = zeros(1,n_uni);
    get_ui = zeros(length(os),max_pi);
    ui = 0;
    for oi = 1:length(os)
        for pi = 1:os(oi).n
            ui = ui+1;
            get_oi(ui) = oi;
            get_pi(ui) = pi;
            get_ui(oi,pi) = ui;
        end
    end
end


%%% read input file and create corresponding origami objects
function [os,linked,o_types,conn_types,conn_r12_eqs] = read_input(inFile)

    %%% open file
    f = fopen(inFile, 'r');
    if f == -1
        error("Could not open file.");
    end

    %%% initialize dictionaries
    blocks = dictionary();
    origamis = dictionary();
    o_blocks = dictionary();
    o_conns = dictionary();
    o_counts = dictionary();

    %%% initialize lists and counters
    os = origami.empty;
    linker_otypes = strings(2,0);
    linker_bis = zeros(2,0);
    linker_locs = cell(2,0);
    nlinker = 0;

    %%% loop over lines, discarding empty ones and comments
    linker_flag = false;
    while ~feof(f)
        line = fgetl(f);
        if isempty(line) || startsWith(line, '%')
            continue;
        end

        %%% extract values before comment
        extract = split(extractBefore([line,'%'],'%'));
        extract = extract(~cellfun('isempty',extract));
        if length(extract) < 2
            continue;
        end

        %%% set the appropriate value
        if linker_flag == true
            linker_flag = false;
            linker_otypes(2,nlinker) = extract{1};
            if extract{2} == 'A'
                linker_bis(2,nlinker) = 0;
            elseif extract{2} == 'B'
                linker_bis(2,nlinker) = -1;
            elseif ~isnan(str2double(extract{2}))
                linker_bis(2,nlinker) = str2double(extract{2});
            else
                error("Unknown linker block: " + extract{2} + ".")
            end
            linker_locs{2,nlinker} = extract{3};
        else
            switch extract{1}
                case 'block'
                    blocks(extract{2}) = block(convertCharsToStrings(extract{3}),str2double(extract{4}));
                case 'patch'
                    blocks(extract{3}) = blocks(extract{3}).add_patch(extract{2},str2double(extract{4}),str2double(extract{5}),str2double(extract{6}));
                case 'origami'
                    origamis(extract{2}) = origami();
                    o_blocks{extract{2}} = block.empty;
                    o_conns{extract{2}} = cell(0);
                    o_counts(extract{2}) = 0;
                case 'linker'
                    linker_flag = true;
                    nlinker = nlinker + 1;
                    linker_otypes(1,nlinker) = extract{2};
                    if extract{3} == 'A'
                        linker_bis(1,nlinker) = 0;
                    elseif extract{3} == 'B'
                        linker_bis(1,nlinker) = -1;
                    elseif ~isnan(str2double(extract{3}))
                        linker_bis(1,nlinker) = str2double(extract{3});
                    else
                        error("Unknown linker block: " + extract{3} + ".")
                    end
                    linker_locs{1,nlinker} = extract{4};
                otherwise
                    if isKey(o_counts,extract{1})
                        switch extract{2}
                            case 'blocks'
                                block_list = block.empty;
                                for bi = 1:length(extract)-2
                                    block_list(bi) = blocks(extract{bi+2});
                                end
                                o_blocks{extract{1}} = block_list;
                            case 'conn'
                                nconn = length(o_conns{extract{1}});
                                o_conns{extract{1}}{nconn+1} = {str2double(extract{3}),extract{4},str2double(extract{5}),extract{6},str2double(extract{7})};
                            case 'count'
                                o_counts(extract{1}) = str2double(extract{3});
                            otherwise
                                error("Unknown origami parameter: " + extract{2})
                        end
                    else
                        error("Unknown system parameter: " + extract{1})
                    end
            end
        end
    end
    fclose(f);

    %%% create type array
    o_counts_values = o_counts.values();
    o_types = zeros(1,sum(o_counts_values));
    o_type = 1;
    for oi = 1:length(o_types)
        if oi > sum(o_counts_values(1:o_type))
            o_type = o_type + 1;
        end
        o_types(oi) = o_type;
    end

    %%% create origamis
    o_indices = dictionary();
    o_keys = keys(o_counts)';
    for o_name = o_keys
        o_indices{o_name} = length(os)+1:length(os)+o_counts(o_name);
        os(o_indices{o_name}) = origami( o_blocks{o_name}, o_conns{o_name});
    end

    %%% create linkers
    linked = zeros(length(os),max([os.n]));
    for li = 1:nlinker
        %%% starting side of the linker
        oi1s = o_indices{linker_otypes(1,li)};
        if linker_bis(1,li) == 0
            bi1s = 1:length([os(oi1s(1)).bs]);
        elseif linker_bis(1,li) == -1
            bi1s = 2:length([os(oi1s(1)).bs]);
        else
            bi1s = linker_bis(1,li);
        end
        loc1 = linker_locs{1,li};
        %%% ending side of the linker
        oi2s = o_indices{linker_otypes(2,li)};
        if linker_bis(2,li) == 0
            bi2s = 1:length([os(oi2s(1)).bs]);
        elseif linker_bis(2,li) == -1
            bi2s = 2:length([os(oi2s(1)).bs]);
        else
            bi2s = linker_bis(2,li);
        end
        loc2 = linker_locs{2,li};
        %%% create linker
        linked = linker(oi1s,bi1s,loc1,oi2s,bi2s,loc2,os,linked);
    end

    %%% connections
    no_type = numEntries(o_counts);
    max_nconn = 0;
    for o_name = o_keys
        nconn = length(o_conns{o_name});
        max_nconn = max([nconn,max_nconn]);
    end
    conn_types = zeros(no_type,max_nconn);
    conn_r12_eqs = zeros(no_type,max_nconn);
    conn_type_count = 0;
    o_type = 0;
    for o_name = o_keys
        o_type = o_type+1;
        for j = 1:length(o_conns{o_name})
            conn_type_count = conn_type_count+1;
            conn_types(o_type,j) = conn_type_count;
            conn_r12_eqs(o_type,j) = o_conns{o_name}{j}{5};
        end
    end
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% File Writers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% write lammps geometry file
function compose_geo(geoFile,geoVisFile,os,linked,o_types,conn_types,dbox)

    %%% grouping atoms
    % atom ID - universal bead index
    % molecule tag - universal block index
    % atom type - unlinked bead (1), unlinked patch (2), linker bead/patch (>2)
    % bond type - linker (1), connection (>1)

    %%% calculate index mapping
    [get_oi,get_pi,get_ui] = map_indices(os);
    nlinker = ars.myMax(linked)/2;
    
    %%% unwrap from pbc
    for oi = 1:length(os)
        for bi = 1:length(os(oi).bs)
            r_ref = ars.calcCOM(os(oi).bs(bi).r,dbox);
            for i = 1:os(oi).bs(bi).n
                pi = os(oi).get_pi(bi,i);
                os(oi).r(:,pi) = r_ref + ars.applyPBC(os(oi).bs(bi).r(:,i) - r_ref, dbox);
            end
        end
    end

    %%% get count for objects
    natom = sum([os.n]);
    nbond = 0;
    for oi = 1:length(os)
        nbond = nbond + length(os(oi).cs);
    end

    %%% compile atom info
    atoms = zeros(5,natom);
    linker_isPatch = zeros(nlinker);
    block_count = 0;
    for ui = 1:natom
        oi = get_oi(ui);
        pi = get_pi(ui);
        if os(oi).get_i(pi) == 1
            block_count = block_count + 1;
        end
        atoms(1,ui) = block_count;
        if linked(oi,pi) > 0
            atoms(2,ui) = linked(oi,pi) + 2;
            linker_index = ceil(linked(oi,pi)/2);
            if linker_isPatch(linker_index) == 0
                if ~os(oi).is_real(pi)
                    linker_isPatch(linker_index) = 1;
                end
            end
        elseif os(oi).get_i(pi) > os(oi).bs(os(oi).get_bi(pi)).n_r
            atoms(2,ui) = 2;
        else
            atoms(2,ui) = 1;
        end
        atoms(3:5,ui) = os(oi).r(:,pi);
    end

    %%% compile bond info
    bonds = zeros(3,nbond);
    if nbond > 0
        bond_count = 0;
        for oi = 1:length(os)
            for ci = 1:length(os(oi).cs)
                bond_count = bond_count + 1;
                ui1 = get_ui( oi, os(oi).get_pi( os(oi).cs(ci).bis(1), os(oi).cs(ci).is(1) ) );
                ui2 = get_ui( oi, os(oi).get_pi( os(oi).cs(ci).bis(2), os(oi).cs(ci).is(2) ) );
                bond_type = conn_types(o_types(oi),ci) + 1;
                bonds(1,bond_count) = bond_type;
                bonds(2,bond_count) = ui1;
                bonds(3,bond_count) = ui2;
            end
        end
    end

    %%% compile angle info
    angles = zeros(4,0);

    %%% compile mass info
    natomType = ars.myMax(linked) + 2;
    masses = ones(1,natomType);
    masses(2) = 0.01;
    for i = 3:natomType
        linker_index = ceil((i-2)/2);
        if linker_isPatch(linker_index) == 1
            masses(i) = 0.01;
        end
    end

    %%% write simulation geometry file
    ars.writeGeo(geoFile,dbox,atoms,bonds,angles,masses=masses)

    %%% compile atom info for visualization
    atoms_vis = atoms;
    no_type = max(o_types);
    for ui = 1:natom
        oi = get_oi(ui);
        pi = get_pi(ui);
        atoms_vis(1,ui) = oi;
        if linked(oi,pi) > 0
            atoms_vis(2,ui) = no_type + 1 + floor((linked(oi,pi)+1)/2);
        elseif os(oi).get_i(pi) > os(oi).bs(os(oi).get_bi(pi)).n_r
            atoms_vis(2,ui) = no_type + 1;
        else
            atoms_vis(2,ui) = o_types(oi);
        end
    end

    %%% write visualization geometry file
    ars.writeGeo(geoVisFile,dbox,atoms_vis,bonds,angles)
end


%%% write lammps input file
function write_input(inputFile,p,linked,conn_types,conn_r12_eqs)

    %%% formatting
    nlinker = ars.myMax(linked)/2;
    natomType = 2+nlinker*2;
    len_nlinker = floor(log10(nlinker))+1;
    len_natomType = floor(log10(natomType))+1;

    %%% calculations
    comm_cutoff = max([p.r12_cut_WCA, p.r12_eq_linker]);

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
    fprintf(f,strcat(...
        "bond_style      harmonic\n",...
        "bond_coeff      1 ", ars.fstring(p.k_x_linker/2,0,2), " ", ars.fstring(p.r12_eq_linker,0,2), "\n"));
    no_type = size(conn_types,1);
    for i = 1:no_type
        nconn_type = nnz(conn_types(i,:));
        for j = 1:nconn_type
            fprintf(f,strcat(...
                "bond_coeff      ", num2str(conn_types(i,j)+1), " ", ars.fstring(p.k_x_conn/2,0,2), " ", ars.fstring(conn_r12_eqs(i,j),0,2), "\n"));
        end
    end
    
    %%% define groups
    if nlinker > 0
        fprintf(f,strcat(...
            "group           linked type 3:", num2str(natomType), "\n"));
    end
    fprintf(f,"\n");

    %%% setup simulation
    fprintf(f,strcat(...
        "## Thermostat\n",...
        "fix             tstat all rigid/nve molecule langevin ", num2str(p.T), " ", num2str(p.T), " 0.049 37\n",...
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

    %%% production settings
    if p.nstep_prod > 0
        fprintf(f,strcat(...
            "## Updates\n",...
            "dump            dump1 all custom ", num2str(p.dump_every), " trajectory.dat id mol xs ys zs\n",...
            "dump_modify     dump1 sort id\n",...
            "compute         comp1 all property/local btype batom1 batom2\n",...
            "dump            dump2 all local ", num2str(p.dump_every), " dump_bonds.txt index c_comp1[1] c_comp1[2] c_comp1[3]\n\n"));
    end

    %%% define linkers
    if nlinker > 0
        fprintf(f,"## Reactions\n");
        for li = 1:nlinker
            fixName = strcat("linker",ars.fstring(li,len_nlinker,0,"L"));
            add_bond_create(f,fixName,li,p.r12_eq_linker,p.react_every,len_natomType);
        end
        fprintf(f,"\n");
    end

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
function add_bond_create(f,fixName,li,r12_cut,react_every,len_natomType)
    fprintf(f,strcat(...
        "fix             ", fixName, " linked bond/create ", num2str(react_every), " ",...
        ars.fstring(li*2+1,len_natomType,0,"L"), " ", ars.fstring(li*2+2,len_natomType,0,"L"), " ",...
        ars.fstring(r12_cut,0,2), " 1 ",...
        "iparam 1 ", ars.fstring(li*2+1,len_natomType,0,"L"), " ",...
        "jparam 1 ", ars.fstring(li*2+2,len_natomType,0,"L"), "\n"));
end
