%%% Housekeeping
clc; clear; close all;

%%% Description
% this script simulates and/or writes the LAMMPS scripts necessary for
  % simulating block/tether-style DNA origami (rigid bodies connected by
  % flexible chains, with selective sticky patches).

%%% Version Description
% added capabilities: better input file

%%% To Do
% remove tether? if not, add capability for connection to handle tethers

%%% Notation
% r - position vector
% r12 - vector pointing from position of bead 1 to bead 2
% r_h - hydrodynamic radius
% o, os, oi - origami, list of origamis, origami index
% t, ts, ti - tether, list of tethers, tether index
% b, bs, bi - block, list of blocks, block index
% nvar - number of var
% var.n or n_var - number of beads in var (e.g. t.n = beads in tether)
% i - bead index within domains (tethers, blocks)
% pi - bead (particle) index across entire origami (tethers then blocks)
% pai - patch index (within block)
% connection - permenant (usually ssDNA scaffold) bond
% linker - switchable (usually hybridizing DNA) bond

%%% Input file
% blocks: rigid bodies
  % name - how to identify the block for adding patches and making origamis
  % pattern - arrangement of beads in the xy-plane (see block object)
  % height - number of beads in the z-direction
% tethers: womrm-like-chains
  % name - how to identify the tether for adding patches and making origamis
  % length - number of beads in the chain
% patches: locations on blocks used for linkers
  % block type - name of block on which to place the patch
  % theta - polar angle in xy-plane of patch location
  % radius - distance (in r12_adj_block units) from z-axis of patch location
  % z - height (in r12_eq_block units) along z-axis of patch location
% origami: connected collection of blocks and tethers
  % name - how to identify the origami for adding linkers
  % design - arrangement of blocks and tethers (see origami object)
  % block - names of block types to use
  % conn - location for 5' and 3' connections
  % count - number of origamis to create
% linker: bonds that can form if beads get close enough
  % origami type - name of origami on which to create the linker
  % block indices - number, "A" for all blocks, or "B" for all but the last block
  % conn - location for 5' and 3' connections
% note: locations are defined by a flag (P or B, for patch or bead),
  % followed by a hyphen, followed by either the patch index or the bead
  % location, as given by the helix identifier (L, M, or R, for left, 
  % middle, or right) and the z-location; this applies to locations given
  % for origami connections and linkers.
% note: linkers require two lines, one starting with the keywork "linker"
  % that defines the location of the first side of the linker, and another
  % directly afterwards that defines the location of the second side of the
  % linker.


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% read input
inFile = "./designs/dimer3C.txt";
[os,linked,is_intra,o_types,conn_types,conn_r12_eqs] = read_input(inFile);

%%% output parameters
simsFold = "/Users/dduke/Files/block_tether/free_origami/lammps/active/";
outputFold = "";
nsim = 1;

%%% simulation parameters
nstep_eq        = 1E4;      %steps    - if and how long to equilibrate/shrink
shrink_ratio    = 1;        %none     - box compression (final/initial)
nstep_prod      = 1E7;      %steps    - if and how long to run producton
dump_every      = 1E4;      %steps    - how often to write to output

%%% computational parameters
dt              = 0.08;     %ns       - timestep
dbox            = 150;      %nm       - periodic boundary diameter
verlet_skin     = 4;        %nm       - width of neighbor list skin (= r12_cut - r12_cut_WCA)
neigh_every     = 1E1;      %steps    - how often to consider updating neighbor list
react_every     = 1E1;      %steps    - how often to check for linker hybridization

%%% physical parameters
kB              = 0.0138;   %pN*nm/K  - Boltzmann constant
T               = 300;      %K        - temperature
r_h_bead        = 1.28;     %nm       - hydrodynamic radius of single bead
visc            = 0.8472;   %mPa/s    - viscosity (pN*ns/nm^2)

%%% bead parameters
sigma           = 10;       %nm       - bead van der Walls radius
epsilon         = 1;        %pN*nm    - WCA energy parameter
r12_eq_block    = 5;        %nm       - equilibrium block bead separation
r12_adj_block   = 5;        %nm       - adjacent block helix separation

%%% tether parameters 
r12_eq_tether   = 2.72;     %nm       - equilibrium bead separation
k_x_tether      = 152;      %pN/nm    - backbone spring constant
Lp_tether       = 4;        %nm       - persistence length

%%% ghost tether parameters
k_x_ghost       = 5.18;     %pN/nm    - attractive backbone spring constant (5.18 for ds14)

%%% linker parameters
nt_linker       = 14;       %nt       - number of nucleotides
k_x_linker      = 55.68;    %pN/nm    - attractive backbone spring constant (55.68 for ds14)
U_cut_linker    = 10;       %kcal/mol - energy cutoff for linker attraction

%%% create parameters class
p = parameters(nstep_eq,shrink_ratio,nstep_prod,dump_every,...
               dt,dbox,verlet_skin,neigh_every,react_every,kB,T,r_h_bead,visc,...
               sigma,epsilon,r12_eq_block,r12_adj_block,r12_eq_tether,k_x_tether,Lp_tether,...
               k_x_ghost,nt_linker,k_x_linker,U_cut_linker);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Write LAMMPS Files %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% output folder 
mkdir(simsFold)

for i = 1:nsim
    if nsim == 1
        outFold = simsFold;
    else
        outFold = simsFold + "sim" + ars.fstring(i,2,0,"R","zero") + "/";
        mkdir(outFold)
    end

    %%% initialize positions
    max_attempts = 10;
    os = place_origami(os,p,max_attempts);

    %%% write lammps simulation geometry file
    geoFile = outFold + "geometry.in";
    geoVisFile = outFold + "geometry_vis.in";
    [nbond,nangle] = compose_geo(geoFile,geoVisFile,os,linked,o_types,conn_types,dbox);
    
    %%% write lammps input file
    inputFile = outFold + "lammps.in";
    write_input(inputFile,outputFold,p,is_intra,nbond,nangle,conn_r12_eqs)
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Simulation Functions %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% initilize positions of the entire system
function os = place_origami(os,p,max_attempts)
    for attempts = 1:max_attempts
        r_other = [];
        for oi = 1:length(os)
            [os(oi),overlap] = os(oi).init(p,r_other);
            if overlap == true
                if attempts >= max_attempts
                    error("Could not place origamis.")
                end
                continue
            end
            os(oi).r = os(oi).update_r;
            r_other = ars.my_horzcat(r_other,os(oi).r);
        end
        if overlap == false
            break
        end
    end
    fprintf("Initialization complete.\n")
end


%%% create linkers by updating linked matrix
function [linked,is_intra] = linker(oi1s,bi1s,loc1,oi2s,bi2s,loc2,os,linked,is_intra)

    %%% meaning of linked matrix values
    % 0 - not linked to anything
    % odd - linked to li+1
    % even - linked to li-1

    %%% identify if linker(s) are intra or inter
    if isequal(oi1s,oi2s)
        interlink = false;
        nlinker = length(is_intra);
        is_intra = [is_intra,ones(1,length(oi1s))];
    else
        interlink = true;
        nlinker = length(is_intra) + 1;
        is_intra = [is_intra,0];
    end

    %%% set values in linker array
    for oi1 = oi1s
        for oi2 = oi2s
            if oi2 == oi1 || interlink == true
                if interlink == false
                    nlinker = nlinker + 1;
                end
                for j1 = 1:length(bi1s)
                    for j2 = 1:length(bi2s)
                        bi1 = bi1s(j1);
                        bi2 = bi2s(j2);
                        i1 = os(oi1).bs(bi1).interpret_loc(loc1);
                        i2 = os(oi2).bs(bi2).interpret_loc(loc2);
                        pi1 = os(oi1).get_pi(2,bi1,i1);
                        pi2 = os(oi2).get_pi(2,bi2,i2);
                        linked(oi1,pi1) = nlinker*2-1;
                        linked(oi2,pi2) = nlinker*2;
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
function [os,linked,is_intra,o_types,conn_types,conn_r12_eqs] = read_input(inFile)

    %%% open file
    f = fopen(inFile, 'r');
    if f == -1
        error("Could not open file");
    end

    %%% initialize dictionaries
    tethers = dictionary();
    blocks = dictionary();
    o_designs = dictionary();
    o_tethers = dictionary();
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
            if extract{2} == "A"
                linker_bis(2,nlinker) = 0;
            elseif extract{2} == "B"
                linker_bis(2,nlinker) = -1;
            elseif ~isnan(str2double(extract{2}))
                linker_bis(2,nlinker) = str2double(extract{2});
            else
                error("Unknown linker block: " + extract{2})
            end
            linker_locs{2,nlinker} = extract{3};
        else
            switch extract{1}
                case 'tether'
                    tethers(extract{2}) = tether(str2double(extract{3}));
                case 'block'
                    blocks(extract{2}) = block(convertCharsToStrings(extract{3}),str2double(extract{4}));
                case 'patch'
                    blocks(extract{3}) = blocks(extract{3}).add_patch(extract{2},str2double(extract{4}),str2double(extract{5}),str2double(extract{6}));
                case 'origami'
                    o_designs(extract{2}) = "";
                    o_tethers{extract{2}} = tether.empty;
                    o_blocks{extract{2}} = block.empty;
                    o_conns{extract{2}} = cell(0);
                    o_counts(extract{2}) = 0;
                case 'linker'
                    linker_flag = true;
                    nlinker = nlinker + 1;
                    linker_otypes(1,nlinker) = extract{2};
                    if extract{3} == "A"
                        linker_bis(1,nlinker) = 0;
                    elseif extract{3} == "B"
                        linker_bis(1,nlinker) = -1;
                    elseif ~isnan(str2double(extract{3}))
                        linker_bis(1,nlinker) = str2double(extract{3});
                    else
                        error("Unknown linker block: " + extract{3})
                    end
                    linker_locs{1,nlinker} = extract{4};
                otherwise
                    if isKey(o_designs,extract{1})
                        switch extract{2}
                            case 'design'
                                o_designs(extract{1}) = extract{3};
                            case 'tether'
                                tether_list = tether.empty;
                                for ti = 1:length(extract)-2
                                    tether_list(ti) = tethers(extract{ti+2});
                                end
                                o_tethers{extract{1}} = tether_list;
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
                                warning("Unknown origami parameter: " + extract{2})
                        end
                    else
                        warning("Unknown system parameter: " + extract{1})
                    end
            end
        end
    end
    fclose(f);

    %%% create type array
    o_counts_values = o_counts.values();
    o_types = zeros(1,sum(o_counts_values));
    type_index = 1;
    for oi = 1:length(o_types)
        if oi > sum(o_counts_values(1:type_index))
            type_index = type_index + 1;
        end
        o_types(oi) = type_index;
    end

    %%% create origamis
    o_indices = dictionary();
    o_keys = keys(o_designs)';
    for o_design = o_keys
        o_indices{o_design} = length(os)+1:length(os)+o_counts(o_design);
        os(o_indices{o_design}) = origami( o_designs(o_design), o_tethers{o_design}, o_blocks{o_design}, o_conns{o_design});
    end

    %%% create linkers
    linked = zeros(length(os),max([os.n]));
    is_intra = zeros(0);
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
        [linked,is_intra] = ...
        linker(oi1s,bi1s,loc1, ...
               oi2s,bi2s,loc2,os,linked,is_intra);
    end

    %%% connections
    no_type = length(o_designs);
    max_nconn = 0;
    for o_design = o_keys
        nconn = length(o_conns{o_design});
        max_nconn = max([nconn,max_nconn]);
    end
    conn_types = zeros(no_type,max_nconn);
    conn_r12_eqs = zeros(no_type,max_nconn);
    conn_type_count = 0;
    i = 0;
    for o_design = o_keys
        i = i+1;
        for j = 1:length(o_conns{o_design})
            conn_type_count = conn_type_count+1;
            conn_types(i,j) = conn_type_count;
            conn_r12_eqs(i,j) = o_conns{o_design}{j}{5};
        end
    end
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% File Writers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% write lammps geometry file
function [nbond,nangle] = compose_geo(geoFile,geoVisFile,os,linked,o_types,conn_types,dbox)

    %%% grouping atoms
    % atom ID - universal bead index
    % molecule tag - universal domain index
    % atom type - tether (1), unlinked block (2), interlinked block (>2), intralinked block (end)
    % bond type - tether backbone (1), linker (2), ghost tether (>2) 

    %%% calculate index mapping
    [get_oi,get_pi,get_ui] = map_indices(os);
    
    %%% unwrap from pbc
    for oi = 1:length(os)
        for ti = 1:length(os(oi).ts)
            r_ref = ars.calcCOM(os(oi).ts(ti).r,dbox);
            for i = 1:os(oi).ts(ti).n
                pi = os(oi).get_pi(1,ti,i);
                os(oi).r(:,pi) = r_ref + ars.applyPBC(os(oi).ts(ti).r(:,i) - r_ref, dbox);
            end
        end
        for bi = 1:length(os(oi).bs)
            r_ref = ars.calcCOM(os(oi).bs(bi).r,dbox);
            for i = 1:os(oi).bs(bi).n
                pi = os(oi).get_pi(2,bi,i);
                os(oi).r(:,pi) = r_ref + ars.applyPBC(os(oi).bs(bi).r(:,i) - r_ref, dbox);
            end
        end
    end

    %%% get count for objects
    natom = sum([os.n]);
    nbond = 0;
    nangle = 0;
    for oi = 1:length(os)
        n_tether = sum([os(oi).ts.n]);
        nbond = nbond + n_tether - length(os(oi).ts) + length(os(oi).cs);
        nangle = nangle + n_tether - 2*length(os(oi).ts);
    end

    %%% compile atom info
    atoms = zeros(5,natom);
    dom_count = 0;
    for ui = 1:natom
        oi = get_oi(ui);
        pi = get_pi(ui);
        if pi == 1 || os(oi).get_dom(pi) ~= os(oi).get_dom(pi-1) || os(oi).get_domi(pi) ~= os(oi).get_domi(pi-1)
            dom_count = dom_count + 1;
        end
        atoms(1,ui) = dom_count;
        if os(oi).get_dom(pi) == 1
            atoms(2,ui) = 1;
        else
            if linked(oi,pi) > 0
                atoms(2,ui) = linked(oi,pi) + 3;
            elseif os(oi).get_i(pi) > os(oi).bs(os(oi).get_domi(pi)).n - os(oi).bs(os(oi).get_domi(pi)).npatch
                atoms(2,ui) = 3;
            else
                atoms(2,ui) = 2;
            end
        end
        atoms(3:5,ui) = os(oi).r(:,pi);
    end

    %%% compile bond info
    bonds = zeros(3,nbond);
    if nbond > 0
        bond_count = 0;
        for oi = 1:length(os)
            if ~isempty(os(oi).ts)
                for pi = 1:sum([os(oi).ts.n])
                    ui = get_ui(oi,pi);
                    if pi+1 <= sum([os(oi).ts.n])
                        if os(oi).get_domi(pi+1) == os(oi).get_domi(pi)
                            bond_count = bond_count + 1;
                            bonds(1,bond_count) = 1;
                            bonds(2,bond_count) = ui+0;
                            bonds(3,bond_count) = ui+1;
                        end
                    end
                end
            end
            for ci = 1:length(os(oi).cs)
                bond_count = bond_count + 1;
                ui1 = get_ui( oi, os(oi).get_pi( os(oi).cs(ci).doms(1), os(oi).cs(ci).domis(1), os(oi).cs(ci).is(1) ) );
                ui2 = get_ui( oi, os(oi).get_pi( os(oi).cs(ci).doms(2), os(oi).cs(ci).domis(2), os(oi).cs(ci).is(2) ) );
                bond_type = conn_types(o_types(oi),ci) + 2;
                bonds(1,bond_count) = bond_type;
                bonds(2,bond_count) = ui1;
                bonds(3,bond_count) = ui2;
            end
        end
    end

    %%% compile angles info
    angles = zeros(4,nangle);
    if nangle > 0
        angle_count = 0;
        for oi = 1:length(os)
            for pi = 1:sum([os(oi).ts.n])
                ui = get_ui(oi,pi);
                if pi+2 <= sum([os(oi).ts.n])
                    if os(oi).get_domi(pi+2) == os(oi).get_domi(pi)
                        angle_count = angle_count + 1;
                        angles(1,angle_count) = 1;
                        angles(2,angle_count) = ui+0;
                        angles(3,angle_count) = ui+1;
                        angles(4,angle_count) = ui+2;
                    end
                end
            end
        end
    end

    %%% compile mass info
    natomType = ars.my_max(linked) + 3;
    masses = ones(1,natomType);
    for i = 3:natomType
        masses(i) = 0.1;
    end

    %%% write simulation geometry file
    ars.write_geo(geoFile,dbox,atoms,bonds,angles,masses=masses)

    %%% initialize visualization arrays
    atoms_vis = atoms;
    bonds_vis = bonds;
    angles_vis = zeros(4,0);

    %%% compile atom info for visualization
    ntype = max(o_types);
    for ui = 1:natom
        oi = get_oi(ui);
        pi = get_pi(ui);
        atoms_vis(1,ui) = oi;
        if linked(oi,pi) > 0
            atoms_vis(2,ui) = ntype + 1 + floor((linked(oi,pi)+1)/2);
        elseif os(oi).get_i(pi) > os(oi).bs(os(oi).get_domi(pi)).n - os(oi).bs(os(oi).get_domi(pi)).npatch
            atoms_vis(2,ui) = ntype + 1;
        else
            atoms_vis(2,ui) = o_types(oi);
        end
    end

    %%% write visualization geometry file
    ars.write_geo(geoVisFile,dbox,atoms_vis,bonds_vis,angles_vis)
end


%%% write lammps input file
function write_input(inputFile,outputFold,p,is_intra,nbond,nangle,conn_r12_eqs)

    %%% create linker list
    linkers = 1:length(is_intra);

    %%% formatting
    nlinker = length(linkers);
    ntype = nlinker*2+3;
    len_nlinker = floor(log10(nlinker))+1;
    len_ntype = floor(log10(ntype))+1;

    %%% open file
    f = fopen(inputFile,'w');
    
    %%% header
    fprintf(f,strcat(...
        "\n#------ Begin Input ------#\n",...
        "# Written by BTBD.m\n\n"));

    %%% basic setup
    fprintf(f,strcat(...
        "## Initialization\n",...
        "units           nano\n",...
        "dimension       3\n",...
        "boundary        p p p\n",...
        "atom_style      molecular\n",...
        "read_data       geometry.in &\n",...
        "                extra/bond/per/atom 1 &\n",...
        "                extra/special/per/atom 1\n\n"));

    %%% neighbor list
    fprintf(f,strcat(...
        "## System definition\n",...
        "neighbor        ", ars.fstring(p.verlet_skin,0,2), " bin\n",...
        "neigh_modify    every ", num2str(p.neigh_every), "\n"));

    %%% initialize pairwise interactions
    fprintf(f,strcat(...
        "pair_style      hybrid/overlay lj/cut ", ars.fstring(p.r12_cut_WCA,0,2), " zero ", ars.fstring(p.r12_cut_linker,0,2), "\n",...
        "pair_coeff      * * zero\n"));

    %%% excluded volume
    fprintf(f,strcat(...
        "pair_coeff      1*2 1*2 lj/cut ", ars.fstring(p.epsilon,0,2), " ", ars.fstring(p.sigma,0,2), "\n",...
        "pair_modify     pair lj/cut shift yes\n"));

    %%% bonds
    if nbond > 0
    fprintf(f,strcat(...
        "bond_style      harmonic\n",...
        "bond_coeff      1 ", ars.fstring(p.k_x_tether/2,0,2), " ", ars.fstring(p.r12_eq_tether,0,2), "\n",...
        "bond_coeff      2 ", ars.fstring(p.k_x_linker/2,0,2), " ", ars.fstring(p.r12_eq_linker,0,2), "\n"));
    conn_count = 0;
    for i = 1:size(conn_r12_eqs,1)
        for j = 1:size(conn_r12_eqs,2)
            conn_count = conn_count + 1;
            fprintf(f,strcat(...
                "bond_coeff      ", num2str(conn_count+2), " ", ars.fstring(p.k_x_ghost/2,0,2), " ", ars.fstring(conn_r12_eqs(i,j),0,2), "\n"));
        end
    end
    end

    %%% tether angles
    if nangle > 0
    fprintf(f,strcat(...
        "angle_style     harmonic\n",...
        "angle_coeff     1 ", ars.fstring(p.k_theta,0,2), " 180\n"));
    end
    
    %%% define groups
    fprintf(f,strcat(...
        "group           tether type 1\n",...
        "group           block type 2:", num2str(ntype), "\n",...
        "group           patch type 3:", num2str(ntype), "\n",...
        "neigh_modify    exclude molecule/intra block\n\n"));

    %%% setup simulation
    fprintf(f,strcat(...
        "## Simulation Setup\n",...
        "fix             tstat1 block rigid/nve molecule langevin ", num2str(p.T), " ", num2str(p.T), " 0.049 37\n",...
        "fix             tstat2 tether langevin ", num2str(p.T), " ", num2str(p.T), " 0.049 37\n",...
        "fix             tstat3 tether nve\n",...
        "thermo_style    custom step temp pe ke etotal\n",...
        "thermo          ", num2str(p.dump_every), "\n\n"));

    %%% equilibration
    if p.nstep_eq > 0
    fprintf(f,strcat(...
        "## Equilibration Settings\n",...
        "timestep        ", num2str(p.dt/10), "\n"));
    if p.shrink_ratio ~= 1
    fprintf(f,strcat(...
        "fix             deformation all deform 1 &\n",...
        "                x final ", ars.fstring(-p.dbox/2*p.shrink_ratio,0,2), " ", ars.fstring(p.dbox/2*p.shrink_ratio,0,2), " &\n",...
        "                y final ", ars.fstring(-p.dbox/2*p.shrink_ratio,0,2), " ", ars.fstring(p.dbox/2*p.shrink_ratio,0,2), " &\n",...
        "                z final ", ars.fstring(-p.dbox/2*p.shrink_ratio,0,2), " ", ars.fstring(p.dbox/2*p.shrink_ratio,0,2), "\n"));
    end
    fprintf(f,"\n");
    end
    
    %%% run equilibration
    if p.nstep_eq > 0
    fprintf(f,strcat(...
        "## Equilibration Steps\n",...
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
        "## Production Settings\n",...
        "timestep        ", num2str(p.dt), "\n",...
        "dump            dump1 all custom ", num2str(p.dump_every), " ", outputFold, "trajectory.dat id mol xs ys zs\n",...
        "dump_modify     dump1 sort id\n",...
        "compute         comp1 all property/local btype batom1 batom2\n",...
        "dump            dump2 all local ", num2str(p.dump_every), " ", outputFold, "dump_bonds.txt index c_comp1[1] c_comp1[2] c_comp1[3]\n"));
    fprintf(f,"\n");
    end

    %%% run production
    if p.nstep_prod > 0
    fprintf(f,strcat(...
        "## Production Steps\n"));
        for li = linkers
            fixName = strcat("linker",ars.fstring(li,len_nlinker,0,"L"),"_create");
            add_bond_create(f,fixName,li,p.r12_eq_linker,p.react_every,len_ntype);
        end
        if ~isempty(linkers)
        fprintf(f,strcat(...
            "fix             linker_break patch bond/break ", num2str(p.react_every), " 3 ", ars.fstring(p.r12_cut_linker,0,2), "\n"));
        end
        fprintf(f,strcat(...
        "run             ", num2str(p.nstep_prod), "\n"));
        for li = linkers
            fixName = strcat("linker",ars.fstring(li,len_nlinker,0,"L"),"_create");
            fprintf(f,strcat("unfix           ",fixName,"\n"));
        end
        if ~isempty(linkers)
        fprintf(f,strcat(...
            "unfix           linker_break\n"));
        end
    end
    
    %%% finalize and close file
    fprintf(f,"\n#------- End Input -------#\n\n");
    fclose(f);
end


%%% write fix bond create command
function add_bond_create(f,fixName,li,r12_cut,react_every,len_ntype)
    fprintf(f,strcat(...
        "fix             ", fixName, " patch bond/create ", num2str(react_every), " ",...
        ars.fstring(li*2+2,len_ntype,0,"L"), " ", ars.fstring(li*2+3,len_ntype,0,"L"), " ",...
        ars.fstring(r12_cut,0,2), " 2 ",...
        "iparam 1 ", ars.fstring(li*2+2,len_ntype,0,"L"), " ",...
        "jparam 1 ", ars.fstring(li*2+3,len_ntype,0,"L"), "\n"));
end
