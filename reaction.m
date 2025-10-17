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
        nreact          % number of contituent reactions
        n_site5         % number of beads in 5p site
        tis_site5       % atom types of beads in 5p site
        n_site3         % number of beads in 3p site
        tis_site3       % atom types of beads in 3p site
        bonds_init      % pre-reaction bonds
        angles_init     % pre-reaction angles
        nsite5          % number of 5p sites
        site5s_oi       % 5p sites origami index
        site5s_bi       % 5p sites block index
        site5s_ibs      % 5p sites patches bead index
        nsite3          % number of 3p sites
        site3s_oi       % 3p sites origami index
        site3s_bi       % 3p sites block index
        site3s_ibs      % 3p sites patches bead index
    end

    methods
        %%% constructor
        function rxn = reaction(label,style,params,ti_start)
            if nargin > 0
                rxn.label = label;
                rxn.style = style;
                rxn.params = params;
                rxn.ti_start = ti_start;
                [rxn.nti,rxn.nreact,tis_internal_site5,tis_internal_site3] = reaction.init(style);
                [rxn.r12s_max,rxn.r12s_min] = reaction.init_r12s(style,params);
                rxn = rxn.set_tis(tis_internal_site5,tis_internal_site3);
                rxn = rxn.set_geo_init();
                rxn.nsite5 = 0;
                rxn.site5s_oi = [];
                rxn.site5s_bi = [];
                rxn.site5s_ibs = [];
                rxn.nsite3 = 0;
                rxn.site3s_oi = [];
                rxn.site3s_bi = [];
                rxn.site3s_ibs = [];
            end
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%% set atom types
        function rxn = set_tis(rxn,tis_internal_site5,tis_internal_site3)
            rxn.n_site5 = length(tis_internal_site5);
            rxn.tis_site5 = zeros(rxn.n_site5,1);
            for ib = 1:rxn.n_site5
                if tis_internal_site5(ib) == 0
                    rxn.tis_site5(ib) = 2;
                else
                    rxn.tis_site5(ib) = rxn.ti_start-1 + tis_internal_site5(ib);
                end
            end
            rxn.n_site3 = length(tis_internal_site3);
            rxn.tis_site3 = zeros(rxn.n_site3,1);
            for ib = 1:rxn.n_site3
                if tis_internal_site3(ib) == 0
                    rxn.tis_site3(ib) = 2;
                else
                    rxn.tis_site3(ib) = rxn.ti_start-1 + tis_internal_site3(ib);
                end
            end
        end


        %%% initialize bonds that connect all atoms within site
        function rxn = set_geo_init(rxn)
            nbond_site5 = (rxn.n_site5*(rxn.n_site5-1))/2;
            bonds_site5 = ones(nbond_site5,3);
            bond_count = 0;
            for i = 1:rxn.n_site5
                for j = i+1:rxn.n_site5
                    bond_count = bond_count + 1;
                    bonds_site5(bond_count,2) = i;
                    bonds_site5(bond_count,3) = j;
                end
            end
            nbond_site3 = (rxn.n_site3*(rxn.n_site3-1))/2;
            bonds_site3 = ones(nbond_site3,3);
            bond_count = 0;
            for i = 1:rxn.n_site3
                for j = i+1:rxn.n_site3
                    bond_count = bond_count + 1;
                    bonds_site3(bond_count,2) = i+rxn.n_site5;
                    bonds_site3(bond_count,3) = j+rxn.n_site5;
                end
            end
            rxn.bonds_init = [bonds_site5;bonds_site3];
            rxn.angles_init = [];
        end


        %%% add reaction site
        function rxn = add_reactant(rxn,is_5p,o,ois,bi,patches)

            %%% check number of patches
            if is_5p
                if length(patches) ~= rxn.n_site5
                    error("Incorrect number of 5p reactant patches.")
                end
            else
                if length(patches) ~= rxn.n_site3
                    error("Incorrect number of 3p reactant patches.")
                end
            end

            %%% 5p site
            if is_5p

                %%% get counts
                norigami = length(ois);
                rxn.nsite5 = rxn.nsite5 + norigami;

                %%% set origami index
                rxn.site5s_oi = [rxn.site5s_oi, ois];

                %%% set block index
                bis = ones(1,norigami)*bi;
                rxn.site5s_bi = [rxn.site5s_bi, bis];
    
                %%% set bead index for each patch
                ibs = zeros(rxn.n_site5,norigami);
                for is = 1:rxn.n_site5
                    ibs(is,:) = o.bs(bi).get_ib_patch(patches(is));
                end
                rxn.site5s_ibs = [rxn.site5s_ibs, ibs];
            
            %%% 3p site
            else

                %%% get counts
                norigami = length(ois);
                rxn.nsite3 = rxn.nsite3 + norigami;

                %%% set origami index
                rxn.site3s_oi = [rxn.site3s_oi, ois];

                %%% set block index
                bis = ones(1,norigami)*bi;
                rxn.site3s_bi = [rxn.site3s_bi, bis];
    
                %%% set bead index for each patch
                ibs = zeros(rxn.n_site3,norigami);
                for is = 1:rxn.n_site3
                    ibs(is,:) = o.bs(bi).get_ib_patch(patches(is));
                end
                rxn.site3s_ibs = [rxn.site3s_ibs, ibs];
            end
        end


        %%% write molecule template file
        function write_molecule(rxn,reactFold,ri,suffix,bonds,angles)

            %%% get counts
            natom = rxn.n_site5 + rxn.n_site3;
            nbond = size(bonds,1);
            nangle = size(angles,1);

            %%% open file
            molFile = reactFold + string(rxn.label) + "_r" + string(ri) + "_" + suffix + ".txt";
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
            
            %%% atoms
            atom_count = 0;
            fprintf(f,strcat(...
                "\nTypes\n\n"));
            for ir = 1:rxn.n_site5
                atom_count = atom_count + 1;
            fprintf(f,strcat(...
                string(atom_count),"\t",string(rxn.tis_site5(ir)),"\n"));
            end
            for ir = 1:rxn.n_site3
                atom_count = atom_count + 1;
            fprintf(f,strcat(...
                string(atom_count),"\t",string(rxn.tis_site3(ir)),"\n"));
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
        end


        %%% write map file
        function write_map(rxn,reactFold,ri,initiators,acons)

            arguments
                rxn
                reactFold
                ri
                initiators
                acons = []
            end
           
            %%% get counts
            natom = rxn.n_site5 + rxn.n_site3;

            %%% count constraints
            ncons = size(acons,1);

            %%% open file
            mapFile = reactFold + string(rxn.label) + "_r" + string(ri) + "_map.txt";
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

            %%% constraints
            if ncons > 0
            fprintf(f,strcat(...
                "\nConstraints\n\n"));
                for aci = 1:size(acons,1)
                fprintf(f,strcat(...
                    "angle ",...
                    string(acons(aci,1))," ",...
                    string(acons(aci,2))," ",...
                    string(acons(aci,3))," ",...
                    string(acons(aci,4))," ",...
                    string(acons(aci,5)),"\n"));
                end
            end
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Define Styles %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%% write react files
        function write_react(rxn,reactFold,pots,apots)

            %%% reaction style
            switch rxn.style

                %%% bond
                case "single"
    
                    %%% put names to parameters
                    bond_type = pots(rxn.params{2}).index;
    
                    %%% initiators
                    initiators = [1 2];

                    %%% set post-reaction bonds
                    bonds_pst = [rxn.bonds_init;
                                 bond_type 1 2];
    
                    %%% write files
                    rxn.write_molecule(reactFold,1,"pre",rxn.bonds_init,rxn.angles_init)
                    rxn.write_molecule(reactFold,1,"pst",bonds_pst,rxn.angles_init)
                    write_map(rxn,reactFold,1,initiators)


                %%% breakable bond
                case "single_rev"
    
                    %%% put names to parameters
                    bond_type = pots(rxn.params{3}).index;

                    %%% initiators
                    initiators = [1 2];

                    %%% set post-reaction bonds
                    bonds_pst = [rxn.bonds_init;
                                 bond_type 1 2];
    
                    %%% write files (forward reaction)
                    rxn.write_molecule(reactFold,1,"pre",rxn.bonds_init,rxn.angles_init)
                    rxn.write_molecule(reactFold,1,"pst",bonds_pst,rxn.angles_init)
                    write_map(rxn,reactFold,1,initiators)

                    %%% write files (reverse reaction)
                    rxn.write_molecule(reactFold,2,"pre",bonds_pst,rxn.angles_init)
                    rxn.write_molecule(reactFold,2,"pst",rxn.bonds_init,rxn.angles_init)
                    write_map(rxn,reactFold,2,initiators)


                %%% bond with angles
                case "single_stiff"
    
                    %%% put names to parameters
                    bond_type = pots(rxn.params{2}).index;
                    angle_type = apots(rxn.params{3}).index;
    
                    %%% initiators
                    initiators = [1 3];
    
                    %%% set post-reaction bonds
                    bonds_pst = [rxn.bonds_init;
                                 bond_type 1 3];
    
                    %%% set post-reaction angles
                    angles_pst = [rxn.angles_init;
                                  angle_type 2 1 3;
                                  angle_type 1 3 4];
    
                    %%% write files
                    rxn.write_molecule(reactFold,1,"pre",rxn.bonds_init,rxn.angles_init)
                    rxn.write_molecule(reactFold,1,"pst",bonds_pst,angles_pst)
                    write_map(rxn,reactFold,1,initiators)

                %%% breakable bond with angles
                case "single_stiff_rev"
    
                    %%% put names to parameters
                    bond_type = pots(rxn.params{3}).index;
                    angle_type = apots(rxn.params{4}).index;

                    %%% initiators
                    initiators = [1 3];
    
                    %%% set post-reaction bonds
                    bonds_pst = [rxn.bonds_init;
                                 bond_type 1 3];
    
                    %%% set post-reaction angles
                    angles_pst = [rxn.angles_init;
                                  angle_type 2 1 3;
                                  angle_type 1 3 4];
    
                    %%% write files (forward reaction)
                    rxn.write_molecule(reactFold,1,"pre",rxn.bonds_init,rxn.angles_init)
                    rxn.write_molecule(reactFold,1,"pst",bonds_pst,angles_pst)
                    write_map(rxn,reactFold,1,initiators)

                    %%% write files (reverse reaction)
                    rxn.write_molecule(reactFold,2,"pre",bonds_pst,angles_pst)
                    rxn.write_molecule(reactFold,2,"pst",rxn.bonds_init,rxn.angles_init)
                    write_map(rxn,reactFold,2,initiators)

                %%% bond with angles added later
                case "single_stiffen"
    
                    %%% put names to parameters
                    bond_type = pots(rxn.params{2}).index;
                    angle_type = apots(rxn.params{3}).index;
                    theta_min = str2double(rxn.params{4});

                    %%% initiators
                    initiators = [1 3];
    
                    %%% set post-reaction bonds
                    bonds_pst = [rxn.bonds_init;
                                 bond_type 1 3];
    
                    %%% write files (react 1)
                    rxn.write_molecule(reactFold,1,"pre",rxn.bonds_init,rxn.angles_init)
                    rxn.write_molecule(reactFold,1,"pst",bonds_pst,rxn.angles_init)
                    write_map(rxn,reactFold,1,initiators)

                    %%% set post-reaction angles
                    angles_pst = [rxn.angles_init;
                                  angle_type 2 1 3;
                                  angle_type 1 3 4];

                    %%% set angle constraints
                    acons = [2 1 3 theta_min 180;
                             1 3 4 theta_min 180];

                    %%% write files (react 1)
                    rxn.write_molecule(reactFold,2,"pre",bonds_pst,rxn.angles_init)
                    rxn.write_molecule(reactFold,2,"pst",bonds_pst,angles_pst)
                    write_map(rxn,reactFold,2,initiators,acons)

                %%% breakable bond with angles added later
                case "single_stiffen_rev"
    
                    %%% put names to parameters
                    bond_type = pots(rxn.params{3}).index;
                    angle_type = apots(rxn.params{4}).index;
                    theta_min = str2double(rxn.params{5});

                    %%% initiators
                    initiators = [1 3];
    
                    %%% set post-reaction bonds
                    bonds_pst = [rxn.bonds_init;
                                 bond_type 1 3];
    
                    %%% write files (bond creation)
                    rxn.write_molecule(reactFold,1,"pre",rxn.bonds_init,rxn.angles_init)
                    rxn.write_molecule(reactFold,1,"pst",bonds_pst,rxn.angles_init)
                    write_map(rxn,reactFold,1,initiators)

                    %%% set post-reaction angles
                    angles_pst = [rxn.angles_init;
                                  angle_type 2 1 3;
                                  angle_type 1 3 4];

                    %%% set angle constraints
                    acons = [2 1 3 theta_min 180;
                             1 3 4 theta_min 180];

                    %%% write files (angle creation)
                    rxn.write_molecule(reactFold,2,"pre",bonds_pst,rxn.angles_init)
                    rxn.write_molecule(reactFold,2,"pst",bonds_pst,angles_pst)
                    write_map(rxn,reactFold,2,initiators,acons)

                    %%% write files (angle break)
                    rxn.write_molecule(reactFold,3,"pre",bonds_pst,rxn.angles_init)
                    rxn.write_molecule(reactFold,3,"pst",rxn.bonds_init,rxn.angles_init)
                    write_map(rxn,reactFold,3,initiators)

                %%% two bonds
                case "double"
    
                    %%% put names to parameters
                    bond_type = pots(rxn.params{2}).index;
 
                    %%% set post-reaction bonds
                    bonds_pst = [rxn.bonds_init;
                                 bond_type 1 3;
                                 bond_type 2 4];

                    %%% initiators
                    initiators = [1 3];
    
                    %%% write files (first pair init)
                    rxn.write_molecule(reactFold,1,"pre",rxn.bonds_init,rxn.angles_init)
                    rxn.write_molecule(reactFold,1,"pst",bonds_pst,rxn.angles_init)
                    write_map(rxn,reactFold,1,initiators)

                    %%% initiators
                    initiators = [2 4];
    
                    %%% write files (second pair init)
                    rxn.write_molecule(reactFold,2,"pre",rxn.bonds_init,rxn.angles_init)
                    rxn.write_molecule(reactFold,2,"pst",bonds_pst,rxn.angles_init)
                    write_map(rxn,reactFold,2,initiators)
                
                %%% two bonds with angles added later
                case "double_stiffen"
    
                    %%% put names to parameters
                    bond_type = pots(rxn.params{2}).index;
                    angle_type = apots(rxn.params{3}).index;
                    theta_min = str2double(rxn.params{4});
 
                    %%% set post-reaction bonds
                    bonds_pst = [rxn.bonds_init;
                                 bond_type 1 5;
                                 bond_type 3 7];

                    %%% initiators
                    initiators = [1 5];
    
                    %%% write files (first pair init)
                    rxn.write_molecule(reactFold,1,"pre",rxn.bonds_init,rxn.angles_init)
                    rxn.write_molecule(reactFold,1,"pst",bonds_pst,rxn.angles_init)
                    write_map(rxn,reactFold,1,initiators)

                    %%% initiators
                    initiators = [3 7];
    
                    %%% write files (second pair init)
                    rxn.write_molecule(reactFold,2,"pre",rxn.bonds_init,rxn.angles_init)
                    rxn.write_molecule(reactFold,2,"pst",bonds_pst,rxn.angles_init)
                    write_map(rxn,reactFold,2,initiators)

                    %%% set post-reaction angles
                    angles_pst = [rxn.angles_init;
                                  angle_type 2 1 5;
                                  angle_type 1 5 6;
                                  angle_type 4 3 7;
                                  angle_type 3 7 8];

                    %%% set angle constraints
                    acons = [2 1 5 theta_min 180;
                             1 5 6 theta_min 180;
                             4 3 7 theta_min 180;
                             3 7 8 theta_min 180];

                    %%% write files (angle creation)
                    rxn.write_molecule(reactFold,3,"pre",bonds_pst,rxn.angles_init)
                    rxn.write_molecule(reactFold,3,"pst",bonds_pst,angles_pst)
                    write_map(rxn,reactFold,3,initiators,acons)

                %%% four bonds
                case "quad"
    
                    %%% put names to parameters
                    bond_type = pots(rxn.params{2}).index;
 
                    %%% set post-reaction bonds
                    bonds_pst = [rxn.bonds_init;
                                 bond_type 1 5;
                                 bond_type 2 6;
                                 bond_type 3 7;
                                 bond_type 4 8];

                    %%% initiators
                    initiators = [1 5];
    
                    %%% write files (first pair init)
                    rxn.write_molecule(reactFold,1,"pre",rxn.bonds_init,rxn.angles_init)
                    rxn.write_molecule(reactFold,1,"pst",bonds_pst,rxn.angles_init)
                    write_map(rxn,reactFold,1,initiators)

                    %%% initiators
                    initiators = [2 6];
    
                    %%% write files (second pair init)
                    rxn.write_molecule(reactFold,2,"pre",rxn.bonds_init,rxn.angles_init)
                    rxn.write_molecule(reactFold,2,"pst",bonds_pst,rxn.angles_init)
                    write_map(rxn,reactFold,2,initiators)

                    %%% initiators
                    initiators = [3 7];
    
                    %%% write files (third pair init)
                    rxn.write_molecule(reactFold,3,"pre",rxn.bonds_init,rxn.angles_init)
                    rxn.write_molecule(reactFold,3,"pst",bonds_pst,rxn.angles_init)
                    write_map(rxn,reactFold,3,initiators)

                    %%% initiators
                    initiators = [4 8];
    
                    %%% write files (fourth pair init)
                    rxn.write_molecule(reactFold,4,"pre",rxn.bonds_init,rxn.angles_init)
                    rxn.write_molecule(reactFold,4,"pst",bonds_pst,rxn.angles_init)
                    write_map(rxn,reactFold,4,initiators)

                %%% four bonds with angles added later
                case "quad_stiffen"
    
                    %%% put names to parameters
                    bond_type = pots(rxn.params{2}).index;
                    angle_type = apots(rxn.params{3}).index;
                    theta_min = str2double(rxn.params{4});
 
                    %%% set post-reaction bonds
                    bonds_pst = [rxn.bonds_init;
                                 bond_type 1 9;
                                 bond_type 3 11;
                                 bond_type 5 13;
                                 bond_type 7 15];

                    %%% initiators
                    initiators = [1 9];
    
                    %%% write files (first pair init)
                    rxn.write_molecule(reactFold,1,"pre",rxn.bonds_init,rxn.angles_init)
                    rxn.write_molecule(reactFold,1,"pst",bonds_pst,rxn.angles_init)
                    write_map(rxn,reactFold,1,initiators)

                    %%% initiators
                    initiators = [3 11];
    
                    %%% write files (second pair init)
                    rxn.write_molecule(reactFold,2,"pre",rxn.bonds_init,rxn.angles_init)
                    rxn.write_molecule(reactFold,2,"pst",bonds_pst,rxn.angles_init)
                    write_map(rxn,reactFold,2,initiators)

                    %%% initiators
                    initiators = [5 13];
    
                    %%% write files (third pair init)
                    rxn.write_molecule(reactFold,3,"pre",rxn.bonds_init,rxn.angles_init)
                    rxn.write_molecule(reactFold,3,"pst",bonds_pst,rxn.angles_init)
                    write_map(rxn,reactFold,3,initiators)

                    %%% initiators
                    initiators = [7 15];
    
                    %%% write files (fourth pair init)
                    rxn.write_molecule(reactFold,4,"pre",rxn.bonds_init,rxn.angles_init)
                    rxn.write_molecule(reactFold,4,"pst",bonds_pst,rxn.angles_init)
                    write_map(rxn,reactFold,4,initiators)

                    %%% set post-reaction angles
                    angles_pst = [rxn.angles_init;
                                  angle_type 2 1 9;
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
                    rxn.write_molecule(reactFold,5,"pre",bonds_pst,rxn.angles_init)
                    rxn.write_molecule(reactFold,5,"pst",bonds_pst,angles_pst)
                    write_map(rxn,reactFold,5,initiators,acons)

                %%% four bonds with angles added later
                case "quad_stiffen_help"
    
                    %%% put names to parameters
                    bond_type = pots(rxn.params{2}).index;
                    angle_type_H = apots(rxn.params{3}).index;
                    angle_type = apots(rxn.params{4}).index;
                    theta_min = str2double(rxn.params{5});
 
                    %%% set post-reaction bonds
                    bonds_pst = [rxn.bonds_init;
                                 bond_type 1 9;
                                 bond_type 3 11;
                                 bond_type 5 13;
                                 bond_type 7 15];

                    %%% set post-reaction angles
                    angles_H_pst = [rxn.angles_init;
                                  angle_type_H 2 1 9;
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
                    rxn.write_molecule(reactFold,1,"pre",rxn.bonds_init,rxn.angles_init)
                    rxn.write_molecule(reactFold,1,"pst",bonds_pst,angles_H_pst)
                    write_map(rxn,reactFold,1,initiators)

                    %%% initiators
                    initiators = [3 11];
    
                    %%% write files (second pair init)
                    rxn.write_molecule(reactFold,2,"pre",rxn.bonds_init,rxn.angles_init)
                    rxn.write_molecule(reactFold,2,"pst",bonds_pst,angles_H_pst)
                    write_map(rxn,reactFold,2,initiators)

                    %%% initiators
                    initiators = [5 13];
    
                    %%% write files (third pair init)
                    rxn.write_molecule(reactFold,3,"pre",rxn.bonds_init,rxn.angles_init)
                    rxn.write_molecule(reactFold,3,"pst",bonds_pst,angles_H_pst)
                    write_map(rxn,reactFold,3,initiators)

                    %%% initiators
                    initiators = [7 15];
    
                    %%% write files (fourth pair init)
                    rxn.write_molecule(reactFold,4,"pre",rxn.bonds_init,rxn.angles_init)
                    rxn.write_molecule(reactFold,4,"pst",bonds_pst,angles_H_pst)
                    write_map(rxn,reactFold,4,initiators)

                    %%% set post-reaction angles
                    angles_pst = [rxn.angles_init;
                                  angle_type 2 1 9;
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
                    rxn.write_molecule(reactFold,5,"pre",bonds_pst,rxn.angles_init)
                    rxn.write_molecule(reactFold,5,"pst",bonds_pst,angles_pst)
                    write_map(rxn,reactFold,5,initiators,acons)

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
        function [nti,nreact,tis_site5,tis_site3] = init(style)

            %%% reaction style
            switch style

                %%% bond
                case "single"
                    nti = 2;
                    nreact = 1;
                    tis_site5 = 1;
                    tis_site3 = 2;

                %%% breakable bond
                case "single_rev"
                    nti = 2;
                    nreact = 2;
                    tis_site5 = 1;
                    tis_site3 = 2;

                %%% bond with angles
                case "single_stiff"
                    nti = 2;
                    nreact = 1;
                    tis_site5 = [1;0];
                    tis_site3 = [2;0];

                %%% breakable bond with angles
                case "single_stiff_rev"
                    nti = 2;
                    nreact = 2;
                    tis_site5 = [1;0];
                    tis_site3 = [2;0];

                %%% bond with angles added later
                case "single_stiffen"
                    nti = 2;
                    nreact = 2;
                    tis_site5 = [1;0];
                    tis_site3 = [2;0];

                %%% breakable bond with angles added later
                case "single_stiffen_rev"
                    nti = 2;
                    nreact = 3;
                    tis_site5 = [1;0];
                    tis_site3 = [2;0];

                %%% two bonds
                case "double"
                    nti = 4;
                    nreact = 2;
                    tis_site5 = [1;3];
                    tis_site3 = [2;4];

                %%% two bonds with angles added later
                case "double_stiffen"
                    nti = 4;
                    nreact = 3;
                    tis_site5 = [1;0;3;0];
                    tis_site3 = [2;0;4;0];

                %%% four bonds
                case "quad"
                    nti = 8;
                    nreact = 4;
                    tis_site5 = [1;3;5;7];
                    tis_site3 = [2;4;6;8];

                %%% four bonds with angles added later
                case "quad_stiffen"
                    nti = 8;
                    nreact = 5;
                    tis_site5 = [1;0;3;0;5;0;7;0];
                    tis_site3 = [2;0;4;0;6;0;8;0];

                %%% four bonds with helping angles, real angles added later
                case "quad_stiffen_help"
                    nti = 8;
                    nreact = 5;
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