%%% bond potential class for BTBD
classdef reaction
    properties
        label           % reaction name
        style           % reaction style
        params          % other parameters
        ti_start        % starting atom type
        r12s_max        % minimum reaction distances
        r12s_min        % minimum reaction distances
        nti             % number of internal atom types
        nsubreact       % number of contituent reactions
        n_site5         % number of beads in 5p site
        tis_site5       % atom types of beads in 5p site
        n_site3         % number of beads in 3p site
        tis_site3       % atom types of beads in 3p site
        bonds_init      % pre-reaction bonds
        nsite5          % number of 5p sites
        site5s_oi       % 5p sites origami index
        site5s_bi       % 5p sites block index
        site5s_ibs      % 5p sites patches bead index
        nsite3          % number of 3p sites
        site3s_oi       % 3p sites origami index
        site3s_bi       % 3p sites block index
        site3s_ibs      % 3p sites patches bead index
        is_charged      % whether reaction updates charges
    end

    methods
        %%% constructor
        function r = reaction(label,style,params,ti_start)
            if nargin > 0
                r.label = label;
                r.style = style;
                r.params = params;
                r.ti_start = ti_start;
                [r.nti,r.nsubreact,tis_internal_site5,tis_internal_site3,r.is_charged] = reaction.init(style);
                [r.r12s_max,r.r12s_min] = reaction.init_r12s(style,params);
                r = r.set_tis(tis_internal_site5,tis_internal_site3);
                r = r.set_bonds_init();
                r.nsite5 = 0;
                r.site5s_oi = [];
                r.site5s_bi = [];
                r.site5s_ibs = [];
                r.nsite3 = 0;
                r.site3s_oi = [];
                r.site3s_bi = [];
                r.site3s_ibs = [];
            end
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%% set atom types
        function r = set_tis(r,tis_internal_site5,tis_internal_site3)
            r.n_site5 = length(tis_internal_site5);
            r.tis_site5 = zeros(r.n_site5,1);
            for ib = 1:r.n_site5
                if tis_internal_site5(ib) == 0
                    r.tis_site5(ib) = 2;
                else
                    r.tis_site5(ib) = r.ti_start-1 + tis_internal_site5(ib);
                end
            end
            r.n_site3 = length(tis_internal_site3);
            r.tis_site3 = zeros(r.n_site3,1);
            for ib = 1:r.n_site3
                if tis_internal_site3(ib) == 0
                    r.tis_site3(ib) = 2;
                else
                    r.tis_site3(ib) = r.ti_start-1 + tis_internal_site3(ib);
                end
            end
        end


        %%% initialize bonds that connect all atoms within site
        function r = set_bonds_init(r)
            nbond_site5 = (r.n_site5*(r.n_site5-1))/2;
            bonds_site5 = ones(nbond_site5,3);
            bond_count = 0;
            for i = 1:r.n_site5
                for j = i+1:r.n_site5
                    bond_count = bond_count + 1;
                    bonds_site5(bond_count,2) = i;
                    bonds_site5(bond_count,3) = j;
                end
            end
            nbond_site3 = (r.n_site3*(r.n_site3-1))/2;
            bonds_site3 = ones(nbond_site3,3);
            bond_count = 0;
            for i = 1:r.n_site3
                for j = i+1:r.n_site3
                    bond_count = bond_count + 1;
                    bonds_site3(bond_count,2) = i+r.n_site5;
                    bonds_site3(bond_count,3) = j+r.n_site5;
                end
            end
            r.bonds_init = [bonds_site5;bonds_site3];
        end


        %%% add reaction site
        function r = add_site(r,is_5p,o,ois,bi,patches)

            %%% check number of patches
            if is_5p
                if length(patches) ~= r.n_site5
                    error("Incorrect number of patches in 5p site.")
                end
            else
                if length(patches) ~= r.n_site3
                    error("Incorrect number of patches in 3p site.")
                end
            end

            %%% get block indices
            if strcmp(bi,'A')
                bis = 1:length(o.bs);
            elseif strcmp(bi,'B')
                bis = 1:length(o.bs)-1;
            else
                bis = str2double(bi);
            end

            %%% 5p site
            if is_5p
                for oi = ois
                    for bi = bis
                        r.nsite5 = r.nsite5 + 1;
                        r.site5s_oi(r.nsite5) = oi;
                        r.site5s_bi(r.nsite5) = bi;
                        ibs = zeros(r.n_site5,1);
                        for is = 1:r.n_site5
                            ibs(is) = o.bs(bi).get_ib_patch(patches(is));
                        end
                        r.site5s_ibs(:,r.nsite5) = ibs;
                    end
                end
               
            %%% 3p site
            else
                for oi = ois
                    for bi = bis
                        r.nsite3 = r.nsite3 + 1;
                        r.site3s_oi(r.nsite3) = oi;
                        r.site3s_bi(r.nsite3) = bi;
                        ibs = zeros(r.n_site3,1);
                        for is = 1:r.n_site3
                            ibs(is) = o.bs(bi).get_ib_patch(patches(is));
                        end
                        r.site3s_ibs(:,r.nsite3) = ibs;
                    end
                end

            end
        end


        %%% write molecule template file
        function write_mol_file(r,reactFold,ri,is_pre,bonds,angles,dihedrals,options)

            arguments
                r; reactFold; ri; is_pre; bonds
                angles = zeros(0,4)
                dihedrals = zeros(0,5)
                options.charge
                options.initiators
            end

            %%% interpret input
            if r.is_charged
                if ~isfield(options,"charge") || ~isfield(options,"initiators")
                    error("Charge and initiators required to write molecule files for charged reactions.")
                end
                charge = options.charge;
                initiators = options.initiators;
            end

            %%% get counts
            natom = r.n_site5 + r.n_site3;
            nbond = size(bonds,1);
            nangle = size(angles,1);
            ndihedral = size(dihedrals,1);

            %%% open file
            molFile = reactFold + string(r.label) + "_r" + string(ri);
            if is_pre
                molFile = molFile + "_pre.txt";
            else
                molFile = molFile + "_pst.txt";
            end
            f = fopen(molFile,'w');

            %%% header
            fprintf(f,strcat(...
                "## Molecule Template\n",...
                string(natom), " atoms\n"));
            if nbond > 0
            fprintf(f,strcat(...
                string(nbond), " bonds\n"));
            end
            if nangle > 0
            fprintf(f,strcat(...
                string(nangle), " angles\n"));
            end
            if ndihedral > 0
            fprintf(f,strcat(...
                string(ndihedral), " dihedrals\n"));
            end
            if r.is_charged && is_pre
            fprintf(f,strcat(...
                "3 fragments\n"));
            end
            
            %%% atoms
            atom_count = 0;
            fprintf(f,strcat(...
                "\nTypes\n\n"));
            for ir = 1:r.n_site5
                atom_count = atom_count + 1;
            fprintf(f,strcat(...
                string(atom_count),"\t",string(r.tis_site5(ir)),"\n"));
            end
            for ir = 1:r.n_site3
                atom_count = atom_count + 1;
            fprintf(f,strcat(...
                string(atom_count),"\t",string(r.tis_site3(ir)),"\n"));
            end

            %%% bonds
            if nbond > 0
            fprintf(f,strcat(...
                "\nBonds\n\n"));
            end
            for bi = 1:nbond
            fprintf(f,strcat(...
                string(bi),"\t",...
                string(bonds(bi,1)),"\t",...
                string(bonds(bi,2)),"\t",...
                string(bonds(bi,3)),"\n"));
            end

            %%% angles
            if nangle > 0
            fprintf(f,strcat(...
                "\nAngles\n\n"));
            end
            for ai = 1:nangle
            fprintf(f,strcat(...
                string(ai),"\t",...
                string(angles(ai,1)),"\t",...
                string(angles(ai,2)),"\t",...
                string(angles(ai,3)),"\t",...
                string(angles(ai,4)),"\n"));
            end

            %%% dihedrals
            if ndihedral > 0
            fprintf(f,strcat(...
                "\nDihedrals\n\n"));
            end
            for di = 1:ndihedral
            fprintf(f,strcat(...
                string(di),"\t",...
                string(dihedrals(di,1)),"\t",...
                string(dihedrals(di,2)),"\t",...
                string(dihedrals(di,3)),"\t",...
                string(dihedrals(di,4)),"\t",...
                string(dihedrals(di,5)),"\n"));
            end

            %%% charges
            if r.is_charged && ~is_pre
                charges = zeros(natom,1);
                charges(initiators) = charge;
                fprintf(f,strcat(...
                    "\nCharges\n\n"));
                for ir = 1:natom
                fprintf(f,strcat(...
                    string(ir), "\t",...
                    string(charges(ir)), "\n"));
                end
            end

            %%% fragments
            if r.is_charged && is_pre
                fprintf(f,strcat(...
                    "\nFragments\n\n"));
                fprintf(f,strcat(...
                    "1\t", string(initiators(1)), "\n",...
                    "2\t", string(initiators(2)), "\n",...
                    "3\t", string(initiators(1)), " ", string(initiators(2)), "\n"));
            end
        end


        %%% write map file
        function write_map_file(r,reactFold,ri,initiators,acons,dcons,options)

            arguments
                r; reactFold; ri; initiators
                acons double = zeros(0,5)
                dcons double = zeros(0,6)
                options.charge
            end

            %%% interpret input
            if r.is_charged
                if ~isfield(options,"charge")
                    error("Charge required to write map files for charged reactions.")
                end
                charge = options.charge;
            end
           
            %%% get counts
            natom = r.n_site5 + r.n_site3;

            %%% count constraints
            nacons = size(acons,1);
            ndcons = size(dcons,1);
            ncons = nacons + ndcons + r.is_charged*2;

            %%% open file
            mapFile = reactFold + string(r.label) + "_r" + string(ri) + "_map.txt";
            f = fopen(mapFile,'w');
           
            %%% header
            fprintf(f,strcat(...
                "## Reaction Map\n",...
                string(natom), " equivalences\n"));
            if ncons > 0
            fprintf(f,strcat(...
                string(ncons), " constraints\n"));
            end
            fprintf(f,strcat(...
                "\nInitiatorIDs\n\n",...
                string(initiators(1)),"\n",...
                string(initiators(2)),"\n",...
                "\nEquivalences\n\n"));

            %%% atoms
            for ir = 1:natom
            fprintf(f,strcat(...
                string(ir),"\t",string(ir),"\n"));
            end

            %%% constraints header
            if ncons > 0
                fprintf(f,strcat(...
                    "\nConstraints\n\n"));
            end

            %%% angle constraints
            if nacons > 0
                for aci = 1:nacons
                fprintf(f,strcat(...
                    "angle ",...
                    string(acons(aci,1))," ",...
                    string(acons(aci,2))," ",...
                    string(acons(aci,3))," ",...
                    string(acons(aci,4))," ",...
                    string(acons(aci,5)),"\n"));
                end
            end

            %%% dihedral constraints
            if nacons > 0
                for dci = 1:ndcons
                fprintf(f,strcat(...
                    "dihedral ",...
                    string(dcons(dci,1))," ",...
                    string(dcons(dci,2))," ",...
                    string(dcons(dci,3))," ",...
                    string(dcons(dci,4))," ",...
                    string(dcons(dci,5))," ",...
                    string(dcons(dci,6)),"\n"));
                end
            end

            %%% charge constraints
            if r.is_charged
                fprintf(f,strcat(...
                    'custom "round(rxnsum(v_varQ,1)) == ', " ",...
                    string(charge), '"\n'));
                fprintf(f,strcat(...
                    'custom "round(rxnsum(v_varQ,2)) ==', " ",...
                    string(charge), '"\n'));
            end

        end

        
        %%% write lines that define molecule templates
        function write_molecules(r,f)

            %%% loop over sub reactions
            for sri = 1:r.nsubreact

                %%% pre-reaction template
                fprintf(f,strcat(...
                    "molecule        ",...
                    string(r.label), "_r", string(sri), "_pre ",...
                    "react/", string(r.label), "_r" + string(sri) + "_pre.txt\n"));

                %%% post-reaction template
                fprintf(f,strcat(...
                    "molecule        ",...
                    string(r.label), "_r", string(sri), "_pst ",...
                    "react/", string(r.label), "_r" + string(sri) + "_pst.txt\n"));
            end
        end


        %%% write lines that add reactions to fix
        function write_react(r,f,comm_cutoff,react_every)
            for sri = 1:r.nsubreact
                r12_min = r.r12s_min(sri);
                r12_max = r.r12s_max(sri);
                if r12_max == 0
                    r12_max = comm_cutoff;
                end
                fprintf(f,strcat(...
                    " &\n                react ",...
                    string(r.label), "_r", string(sri), " ",...
                    "sticky ", string(react_every), " ",...
                    ars.fstring(r12_min,0,2), " ",...
                    ars.fstring(r12_max,0,2), " ",...
                    string(r.label), "_r", string(sri), "_pre ",...
                    string(r.label), "_r" + string(sri), "_pst ",...
                    "react/", string(r.label), "_r", string(sri), "_map.txt"));
                if r.is_charged
                    fprintf(f," custom_charges 3");
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Define Styles %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%% write react files
        function write_react_files(r,reactFold,pots,apots,dpots)

            %%% reaction style
            switch r.style

                %%% bond
                case "single"
    
                    %%% put names to parameters
                    bond_type = pots(r.params{2}).index;
    
                    %%% initiators
                    initiators = [1 2];

                    %%% set post-reaction bonds
                    bonds_pst = [r.bonds_init;
                                 bond_type 1 2];
    
                    %%% write files
                    r.write_mol_file(reactFold,1,1,r.bonds_init)
                    r.write_mol_file(reactFold,1,0,bonds_pst)
                    write_map_file(r,reactFold,1,initiators)


                %%% breakable bond
                case "single_rev"
    
                    %%% put names to parameters
                    bond_type = pots(r.params{3}).index;

                    %%% initiators
                    initiators = [1 2];

                    %%% set post-reaction bonds
                    bonds_pst = [r.bonds_init;
                                 bond_type 1 2];
    
                    %%% write files (forward reaction)
                    r.write_mol_file(reactFold,1,1,r.bonds_init)
                    r.write_mol_file(reactFold,1,0,bonds_pst)
                    write_map_file(r,reactFold,1,initiators)

                    %%% write files (reverse reaction)
                    r.write_mol_file(reactFold,2,1,bonds_pst)
                    r.write_mol_file(reactFold,2,0,r.bonds_init)
                    write_map_file(r,reactFold,2,initiators)


                %%% bond with angles
                case "single_stiff"
    
                    %%% put names to parameters
                    bond_type = pots(r.params{2}).index;
                    angle_type = apots(r.params{3}).index;
    
                    %%% initiators
                    initiators = [1 3];
    
                    %%% set post-reaction bonds
                    bonds_pst = [r.bonds_init;
                                 bond_type 1 3];
    
                    %%% set post-reaction angles
                    angles_pst = [angle_type 2 1 3;
                                  angle_type 1 3 4];
    
                    %%% write files
                    r.write_mol_file(reactFold,1,1,r.bonds_init)
                    r.write_mol_file(reactFold,1,0,bonds_pst,angles_pst)
                    write_map_file(r,reactFold,1,initiators)

                %%% breakable bond with angles
                case "single_stiff_rev"
    
                    %%% put names to parameters
                    bond_type = pots(r.params{3}).index;
                    angle_type = apots(r.params{4}).index;

                    %%% initiators
                    initiators = [1 3];
    
                    %%% set post-reaction bonds
                    bonds_pst = [r.bonds_init;
                                 bond_type 1 3];
    
                    %%% set post-reaction angles
                    angles_pst = [angle_type 2 1 3;
                                  angle_type 1 3 4];
    
                    %%% write files (forward reaction)
                    r.write_mol_file(reactFold,1,1,r.bonds_init)
                    r.write_mol_file(reactFold,1,0,bonds_pst,angles_pst)
                    write_map_file(r,reactFold,1,initiators)

                    %%% write files (reverse reaction)
                    r.write_mol_file(reactFold,2,1,bonds_pst,angles_pst)
                    r.write_mol_file(reactFold,2,0,r.bonds_init)
                    write_map_file(r,reactFold,2,initiators)

                %%% bond with angles added later
                case "single_stiffen"
    
                    %%% put names to parameters
                    bond_type = pots(r.params{2}).index;
                    angle_type = apots(r.params{3}).index;
                    theta_min = str2double(r.params{4});

                    %%% initiators
                    initiators = [1 3];
    
                    %%% set post-reaction bonds
                    bonds_pst = [r.bonds_init;
                                 bond_type 1 3];
    
                    %%% write files (bond creation)
                    r.write_mol_file(reactFold,1,1,r.bonds_init)
                    r.write_mol_file(reactFold,1,0,bonds_pst)
                    write_map_file(r,reactFold,1,initiators)

                    %%% set post-reaction angles
                    angles_pst = [angle_type 2 1 3;
                                  angle_type 1 3 4];

                    %%% set angle constraints
                    acons = [2 1 3 theta_min 180;
                             1 3 4 theta_min 180];

                    %%% write files (angle creation)
                    r.write_mol_file(reactFold,2,1,bonds_pst)
                    r.write_mol_file(reactFold,2,0,bonds_pst,angles_pst)
                    write_map_file(r,reactFold,2,initiators,acons)

                %%% breakable bond with angles added later
                case "single_stiffen_rev"
    
                    %%% put names to parameters
                    bond_type = pots(r.params{3}).index;
                    angle_type = apots(r.params{4}).index;
                    theta_min = str2double(r.params{5});

                    %%% initiators
                    initiators = [1 3];
    
                    %%% set post-reaction bonds
                    bonds_pst = [r.bonds_init;
                                 bond_type 1 3];
    
                    %%% write files (bond creation)
                    r.write_mol_file(reactFold,1,1,r.bonds_init)
                    r.write_mol_file(reactFold,1,0,bonds_pst)
                    write_map_file(r,reactFold,1,initiators)

                    %%% set post-reaction angles
                    angles_pst = [angle_type 2 1 3;
                                  angle_type 1 3 4];

                    %%% set angle constraints
                    acons = [2 1 3 theta_min 180;
                             1 3 4 theta_min 180];

                    %%% write files (angle creation)
                    r.write_mol_file(reactFold,2,1,bonds_pst)
                    r.write_mol_file(reactFold,2,0,bonds_pst,angles_pst)
                    write_map_file(r,reactFold,2,initiators,acons)

                    %%% write files (angle break)
                    r.write_mol_file(reactFold,3,1,bonds_pst)
                    r.write_mol_file(reactFold,3,0,r.bonds_init)
                    write_map_file(r,reactFold,3,initiators)

                %%% bond with angles and dihedrals added later
                case "single_stiffen_twist"
    
                    %%% put names to parameters
                    bond_type = pots(r.params{2}).index;
                    angle_type = apots(r.params{3}).index;
                    theta_min = str2double(r.params{4});
                    dihedral_type = dpots(r.params{5}).index;
                    phi_max = str2double(r.params{6});

                    %%% initiators
                    initiators = [1 4];
    
                    %%% set post-reaction bonds
                    bonds_pst = [r.bonds_init;
                                 bond_type 1 4];
    
                    %%% define charges
                    charge_pre = 0;
                    charge_pst = 1;
    
                    %%% write files (bond creation)
                    r.write_mol_file(reactFold,1,1,r.bonds_init,charge=charge_pre,initiators=initiators)
                    r.write_mol_file(reactFold,1,0,bonds_pst,charge=charge_pst,initiators=initiators)
                    write_map_file(r,reactFold,1,initiators,charge=charge_pre)

                    %%% set post-reaction angles
                    angles_pst = [angle_type 2 1 4;
                                  angle_type 1 4 5];

                    %%% set angle constraints
                    acons = [2 1 4 theta_min 180;
                             1 4 5 theta_min 180];

                    %%% define charges
                    charge_pre = 1;
                    charge_pst = 2;

                    %%% write files (angle creation)
                    r.write_mol_file(reactFold,2,1,bonds_pst,charge=charge_pre,initiators=initiators)
                    r.write_mol_file(reactFold,2,0,bonds_pst,angles_pst,charge=charge_pst,initiators=initiators)
                    write_map_file(r,reactFold,2,initiators,acons,charge=charge_pre)

                    %%% set post-reaction dihedrals
                    dihedrals_pst = [dihedral_type 3 1 4 6];

                    %%% set angle constraints
                    dcons = [3 1 4 6 -phi_max phi_max];

                    %%% define charges
                    charge_pre = 2;
                    charge_pst = 3;

                    %%% write files (dihedral creation)
                    r.write_mol_file(reactFold,3,1,bonds_pst,angles_pst,charge=charge_pre,initiators=initiators)
                    r.write_mol_file(reactFold,3,0,bonds_pst,angles_pst,dihedrals_pst,charge=charge_pst,initiators=initiators)
                    write_map_file(r,reactFold,3,initiators,acons,dcons,charge=charge_pre)

                %%% two bonds
                case "double"
    
                    %%% put names to parameters
                    bond_type = pots(r.params{2}).index;
 
                    %%% set post-reaction bonds
                    bonds_pst = [r.bonds_init;
                                 bond_type 1 3;
                                 bond_type 2 4];

                    %%% initiators
                    initiators = [1 3];
    
                    %%% write files (first pair init)
                    r.write_mol_file(reactFold,1,1,r.bonds_init)
                    r.write_mol_file(reactFold,1,0,bonds_pst)
                    write_map_file(r,reactFold,1,initiators)

                    %%% initiators
                    initiators = [2 4];
    
                    %%% write files (second pair init)
                    r.write_mol_file(reactFold,2,1,r.bonds_init)
                    r.write_mol_file(reactFold,2,0,bonds_pst)
                    write_map_file(r,reactFold,2,initiators)
                
                %%% two bonds with angles added later
                case "double_stiffen"
    
                    %%% put names to parameters
                    bond_type = pots(r.params{2}).index;
                    angle_type = apots(r.params{3}).index;
                    theta_min = str2double(r.params{4});
 
                    %%% set post-reaction bonds
                    bonds_pst = [r.bonds_init;
                                 bond_type 1 5;
                                 bond_type 3 7];

                    %%% initiators
                    initiators = [1 5];
    
                    %%% write files (first pair init)
                    r.write_mol_file(reactFold,1,1,r.bonds_init)
                    r.write_mol_file(reactFold,1,0,bonds_pst)
                    write_map_file(r,reactFold,1,initiators)

                    %%% initiators
                    initiators = [3 7];
    
                    %%% write files (second pair init)
                    r.write_mol_file(reactFold,2,1,r.bonds_init)
                    r.write_mol_file(reactFold,2,0,bonds_pst)
                    write_map_file(r,reactFold,2,initiators)

                    %%% set post-reaction angles
                    angles_pst = [angle_type 2 1 5;
                                  angle_type 1 5 6;
                                  angle_type 4 3 7;
                                  angle_type 3 7 8];

                    %%% set angle constraints
                    acons = [2 1 5 theta_min 180;
                             1 5 6 theta_min 180;
                             4 3 7 theta_min 180;
                             3 7 8 theta_min 180];

                    %%% write files (angle creation)
                    r.write_mol_file(reactFold,3,1,bonds_pst)
                    r.write_mol_file(reactFold,3,0,bonds_pst,angles_pst)
                    write_map_file(r,reactFold,3,initiators,acons)

                %%% four bonds
                case "quad"
    
                    %%% put names to parameters
                    bond_type = pots(r.params{2}).index;
 
                    %%% set post-reaction bonds
                    bonds_pst = [r.bonds_init;
                                 bond_type 1 5;
                                 bond_type 2 6;
                                 bond_type 3 7;
                                 bond_type 4 8];

                    %%% initiators
                    initiators = [1 5];
    
                    %%% write files (first pair init)
                    r.write_mol_file(reactFold,1,1,r.bonds_init)
                    r.write_mol_file(reactFold,1,0,bonds_pst)
                    write_map_file(r,reactFold,1,initiators)

                    %%% initiators
                    initiators = [2 6];
    
                    %%% write files (second pair init)
                    r.write_mol_file(reactFold,2,1,r.bonds_init)
                    r.write_mol_file(reactFold,2,0,bonds_pst)
                    write_map_file(r,reactFold,2,initiators)

                    %%% initiators
                    initiators = [3 7];
    
                    %%% write files (third pair init)
                    r.write_mol_file(reactFold,3,1,r.bonds_init)
                    r.write_mol_file(reactFold,3,0,bonds_pst)
                    write_map_file(r,reactFold,3,initiators)

                    %%% initiators
                    initiators = [4 8];
    
                    %%% write files (fourth pair init)
                    r.write_mol_file(reactFold,4,1,r.bonds_init)
                    r.write_mol_file(reactFold,4,0,bonds_pst)
                    write_map_file(r,reactFold,4,initiators)

                %%% four bonds with angles added later
                case "quad_stiffen"
    
                    %%% put names to parameters
                    bond_type = pots(r.params{2}).index;
                    angle_type = apots(r.params{3}).index;
                    theta_min = str2double(r.params{4});
 
                    %%% set post-reaction bonds
                    bonds_pst = [r.bonds_init;
                                 bond_type 1 9;
                                 bond_type 3 11;
                                 bond_type 5 13;
                                 bond_type 7 15];

                    %%% initiators
                    initiators = [1 9];
    
                    %%% write files (first pair init)
                    r.write_mol_file(reactFold,1,1,r.bonds_init)
                    r.write_mol_file(reactFold,1,0,bonds_pst)
                    write_map_file(r,reactFold,1,initiators)

                    %%% initiators
                    initiators = [3 11];
    
                    %%% write files (second pair init)
                    r.write_mol_file(reactFold,2,1,r.bonds_init)
                    r.write_mol_file(reactFold,2,0,bonds_pst)
                    write_map_file(r,reactFold,2,initiators)

                    %%% initiators
                    initiators = [5 13];
    
                    %%% write files (third pair init)
                    r.write_mol_file(reactFold,3,1,r.bonds_init)
                    r.write_mol_file(reactFold,3,0,bonds_pst)
                    write_map_file(r,reactFold,3,initiators)

                    %%% initiators
                    initiators = [7 15];
    
                    %%% write files (fourth pair init)
                    r.write_mol_file(reactFold,4,1,r.bonds_init)
                    r.write_mol_file(reactFold,4,0,bonds_pst)
                    write_map_file(r,reactFold,4,initiators)

                    %%% set post-reaction angles
                    angles_pst = [angle_type 2 1 9;
                                  angle_type 1 9 10;
                                  angle_type 4 3 11;
                                  angle_type 3 11 12;
                                  angle_type 6 5 13;
                                  angle_type 5 13 14;
                                  angle_type 8 7 15;
                                  angle_type 7 15 16];

                    %%% set angle constraints
                    acons = [2 1 9   theta_min 180;
                             1 9 10  theta_min 180;
                             4 3 11  theta_min 180;
                             3 11 12 theta_min 180
                             6 5 13  theta_min 180;
                             5 13 14 theta_min 180;
                             8 7 15  theta_min 180;
                             7 15 16 theta_min 180];

                    %%% write files (angle creation)
                    r.write_mol_file(reactFold,5,1,bonds_pst)
                    r.write_mol_file(reactFold,5,0,bonds_pst,angles_pst)
                    write_map_file(r,reactFold,5,initiators,acons)

                %%% four bonds with angles added later
                case "quad_stiffen_help"
    
                    %%% put names to parameters
                    bond_type = pots(r.params{2}).index;
                    angle_type_H = apots(r.params{3}).index;
                    angle_type = apots(r.params{4}).index;
                    theta_min = str2double(r.params{5});
 
                    %%% set post-reaction bonds
                    bonds_pst = [r.bonds_init;
                                 bond_type 1 9;
                                 bond_type 3 11;
                                 bond_type 5 13;
                                 bond_type 7 15];

                    %%% set post-reaction angles
                    angles_H_pst = [angle_type_H 2 1 9;
                                    angle_type_H 1 9 10;
                                    angle_type_H 4 3 11;
                                    angle_type_H 3 11 12;
                                    angle_type_H 6 5 13;
                                    angle_type_H 5 13 14;
                                    angle_type_H 8 7 15;
                                    angle_type_H 7 15 16];

                    %%% initiators
                    initiators = [1 9];
    
                    %%% write files (first pair init)
                    r.write_mol_file(reactFold,1,1,r.bonds_init)
                    r.write_mol_file(reactFold,1,0,bonds_pst,angles_H_pst)
                    write_map_file(r,reactFold,1,initiators)

                    %%% initiators
                    initiators = [3 11];
    
                    %%% write files (second pair init)
                    r.write_mol_file(reactFold,2,1,r.bonds_init)
                    r.write_mol_file(reactFold,2,0,bonds_pst,angles_H_pst)
                    write_map_file(r,reactFold,2,initiators)

                    %%% initiators
                    initiators = [5 13];
    
                    %%% write files (third pair init)
                    r.write_mol_file(reactFold,3,1,r.bonds_init)
                    r.write_mol_file(reactFold,3,0,bonds_pst,angles_H_pst)
                    write_map_file(r,reactFold,3,initiators)

                    %%% initiators
                    initiators = [7 15];
    
                    %%% write files (fourth pair init)
                    r.write_mol_file(reactFold,4,1,r.bonds_init)
                    r.write_mol_file(reactFold,4,0,bonds_pst,angles_H_pst)
                    write_map_file(r,reactFold,4,initiators)

                    %%% set post-reaction angles
                    angles_pst = [angle_type 2 1 9;
                                  angle_type 1 9 10;
                                  angle_type 4 3 11;
                                  angle_type 3 11 12;
                                  angle_type 6 5 13;
                                  angle_type 5 13 14;
                                  angle_type 8 7 15;
                                  angle_type 7 15 16];

                    %%% set angle constraints
                    acons = [2 1 9   theta_min 180;
                             1 9 10  theta_min 180;
                             4 3 11  theta_min 180;
                             3 11 12 theta_min 180
                             6 5 13  theta_min 180;
                             5 13 14 theta_min 180;
                             8 7 15  theta_min 180;
                             7 15 16 theta_min 180];

                    %%% write files (angle stiffening)
                    r.write_mol_file(reactFold,5,1,bonds_pst)
                    r.write_mol_file(reactFold,5,0,bonds_pst,angles_pst)
                    write_map_file(r,reactFold,5,initiators,acons)

                %%% error
                otherwise
                    error("Unknown reaction style.")
            end
        end

    end
    methods (Static)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Static Functions %%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%% initialize reaction parameters based on style
        function [nti,nsubreact,tis_site5,tis_site3,is_charged] = init(style)

            %%% no charge by default
            is_charged = 0;

            %%% reaction style
            switch style

                %%% bond
                case "single"
                    nti = 2;
                    nsubreact = 1;
                    tis_site5 = 1;
                    tis_site3 = 2;

                %%% breakable bond
                case "single_rev"
                    nti = 2;
                    nsubreact = 2;
                    tis_site5 = 1;
                    tis_site3 = 2;

                %%% bond with angles
                case "single_stiff"
                    nti = 2;
                    nsubreact = 1;
                    tis_site5 = [1;0];
                    tis_site3 = [2;0];

                %%% breakable bond with angles
                case "single_stiff_rev"
                    nti = 2;
                    nsubreact = 2;
                    tis_site5 = [1;0];
                    tis_site3 = [2;0];

                %%% bond with angles added later
                case "single_stiffen"
                    nti = 2;
                    nsubreact = 2;
                    tis_site5 = [1;0];
                    tis_site3 = [2;0];

                %%% breakable bond with angles added later
                case "single_stiffen_rev"
                    nti = 2;
                    nsubreact = 3;
                    tis_site5 = [1;0];
                    tis_site3 = [2;0];

                %%% bond with angles and dihedral added later
                case "single_stiffen_twist"
                    nti = 3;
                    nsubreact = 3;
                    tis_site5 = [1;0;3];
                    tis_site3 = [2;0;3];
                    is_charged = 1;

                %%% two bonds
                case "double"
                    nti = 4;
                    nsubreact = 2;
                    tis_site5 = [1;3];
                    tis_site3 = [2;4];

                %%% two bonds with angles added later
                case "double_stiffen"
                    nti = 4;
                    nsubreact = 3;
                    tis_site5 = [1;0;3;0];
                    tis_site3 = [2;0;4;0];

                %%% four bonds
                case "quad"
                    nti = 8;
                    nsubreact = 4;
                    tis_site5 = [1;3;5;7];
                    tis_site3 = [2;4;6;8];

                %%% four bonds with angles added later
                case "quad_stiffen"
                    nti = 8;
                    nsubreact = 5;
                    tis_site5 = [1;0;3;0;5;0;7;0];
                    tis_site3 = [2;0;4;0;6;0;8;0];

                %%% four bonds with helping angles, real angles added later
                case "quad_stiffen_help"
                    nti = 8;
                    nsubreact = 5;
                    tis_site5 = [1;0;3;0;5;0;7;0];
                    tis_site3 = [2;0;4;0;6;0;8;0];
                   
                %%% error
                otherwise
                    error("Unknown reaction type.")
            end
        end

        %%% initialize reaction parameters based on style
        function [r12s_max,r12s_min] = init_r12s(style,params)

            %%% reaction style
            switch style

                %%% bond
                case "single"
                    reaction.check_nparam(params,2)
                    r12s_max = str2double(params{1});
                    r12s_min = 0;

                %%% breakable bond
                case "single_rev"
                    reaction.check_nparam(params,3)
                    r12s_max = [str2double(params{1}) 0];
                    r12s_min = [0 str2double(params{2})];

                %%% bond with angles
                case "single_stiff"
                    reaction.check_nparam(params,3)
                    r12s_max = str2double(params{1});
                    r12s_min = 0;

                %%% breakable bond with angles
                case "single_stiff_rev"
                    reaction.check_nparam(params,4)
                    r12s_max = [str2double(params{1}) 0];
                    r12s_min = [0 str2double(params{2})];

                %%% bond with angles added later
                case "single_stiffen"
                    reaction.check_nparam(params,4)
                    r12_max = str2double(params{1});
                    r12s_max = [r12_max r12_max];
                    r12s_min = [0 0];

                %%% breakable bond with angles added later
                case "single_stiffen_rev"
                    reaction.check_nparam(params,5)
                    r12_max = str2double(params{1});
                    r12s_max = [r12_max r12_max 0];
                    r12s_min = [0 0 str2double(params{2})];

                %%% bond with angles and dihedral added later
                case "single_stiffen_twist"
                    reaction.check_nparam(params,6)
                    r12_max = str2double(params{1});
                    r12s_max = [r12_max r12_max r12_max];
                    r12s_min = [0 0 0];

                %%% two bonds
                case "double"
                    reaction.check_nparam(params,2)
                    r12_max = str2double(params{1});
                    r12s_max = [r12_max r12_max];
                    r12s_min = [0 0];

                %%% two bonds with angles added later
                case "double_stiffen"
                    reaction.check_nparam(params,4)
                    r12_max = str2double(params{1});
                    r12s_max = [r12_max r12_max 0];
                    r12s_min = [0 0 0];

                %%% four bonds
                case "quad"
                    reaction.check_nparam(params,2)
                    r12_max = str2double(params{1});
                    r12s_max = [r12_max r12_max r12_max r12_max];
                    r12s_min = [0 0 0 0];

                %%% four bonds with angles aded later
                case "quad_stiffen"
                    reaction.check_nparam(params,4)
                    r12_max = str2double(params{1});
                    r12s_max = [r12_max r12_max r12_max r12_max 0];
                    r12s_min = [0 0 0 0 0];

                %%% four bonds with angles aded later
                case "quad_stiffen_help"
                    reaction.check_nparam(params,5)
                    r12_max = str2double(params{1});
                    r12s_max = [r12_max r12_max r12_max r12_max 0];
                    r12s_min = [0 0 0 0 0];
                   
                %%% error
                otherwise
                    error("Unknown reaction type.")
            end
        end


        %%% ensure correct number of parameters
        function check_nparam(params,nparam)
            if length(params) ~= nparam
                error("Incorrect number of reaction parameters")
            end
        end

    end
end