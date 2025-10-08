%%% block class for BTBD
classdef block
    properties
        pattern         % helix locations in block xy-plane
        hiH             % helix indices of hollow helices
        hiL             % helix index for left connections
        hiR             % helix index for right connections
        hiM             % helix index for middle connections
        nhelix          % number of helices
        n_helix         % number of beads in each helix
        r_internal      % internal positions (in bead units)
        npatch          % number of patches
        patches         % names of patches
        n_real          % number of real (structural) beads
        n               % number of total beads (including patches)
        r               % bead positions
        get_hi          % map ib to hi
        get_zi          % map ib to zi
        get_ib          % map hi and zi to ib
    end

    methods
        %%% constructor
        function b = block(pattern_label,n_helix)
            if nargin > 0
                [b.pattern,b.hiH,b.hiL,b.hiR,b.hiM] = block.setpat(pattern_label);
                b.nhelix = size(b.pattern,2);
                b.n_helix = n_helix;
                b.npatch = 0;
                b.patches = {};
                b.n_real = b.calc_nbead;
                b.n = b.n_real;
                b.r = zeros(3,b.n);
                [b.get_hi,b.get_zi,b.get_ib] = map_indices(b);
                b.r_internal = init_positions_internal(b);
            end
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%% add patch to block
        function b = add_patch(b,name,x,y,z)
            b.npatch = b.npatch + 1;
            b.patches{b.npatch} = name;
            b.n = b.n + 1;
            b.r(:,b.n) = zeros(3,1);
            b.r_internal(:,b.n) = [x;y;z];
        end


        %%% calcualte number of structural beads in block
        function nbead = calc_nbead(b)
            nbead = size(b.pattern,2)*min(b.n_helix,2) + (size(b.pattern,2)-length(b.hiH))*max(b.n_helix-2,0);
        end


        %%% get block index from geometrical indices
        function [get_hi,get_zi,get_ib] = map_indices(b)
            get_hi = zeros(1,b.n);
            get_zi = zeros(1,b.n);
            get_ib = zeros(b.nhelix,b.n_helix);
            for hi = 1:b.nhelix
                for zi = 1:b.n_helix
                    if zi ~= 1 && zi ~= b.n_helix && any(b.hiH==hi)
                        ib = 0;
                    else
                        ib = size(b.pattern,2)*min(zi-1,1) + ...
                                (size(b.pattern,2)-length(b.hiH))*max(zi-2,0) + hi;
                        if zi ~= 1 && zi ~= b.n_helix
                            ib = ib - sum(b.hiH<hi);
                        end
                    end
                    get_hi(ib) = hi;
                    get_zi(ib) = zi;
                    get_ib(hi,zi) = ib;
                end
            end
        end


        %%% set the positions of the block within its own coordinate system
        function r_internal = init_positions_internal(b)
            r_internal = zeros(3,b.n);
            for ib = 1:b.n_real
                hi = b.get_hi(ib);
                zi = b.get_zi(ib);
                r_internal(1:2,ib) = b.pattern(:,hi);
                r_internal(3,ib) = zi-1;
            end
        end


        %%% get block index from location
        function ib = interpret_loc(b,loc)

            %%% connection to structural bead
            if length(loc) >= 2 && loc(1:2) == "B-"

                %%% find helix
                if loc(3) == 'L'
                    hi = b.hiL;
                elseif loc(3) == 'M'
                    hi = b.hiM;
                elseif loc(3) == 'R'
                    hi = b.hiR;
                else
                    error("Unknown helix type: " + loc(3) + ".")
                end

                %%% find height
                if length(loc) >= 4
                    zi = str2double(loc(4:end));
                else
                    zi = 0;
                end

                %%% get bead index
                ib = b.get_ib(hi,zi);

            %%% connection to patch bead
            else

                %%% find patch
                if loc(1:2) == "P-"
                    patch_index = find(strcmp(b.patches, loc(3:end)));
                else
                    patch_index = find(strcmp(b.patches, loc));
                end

                %%% get bead index
                if ~isempty(patch_index)
                    ib = b.n-b.npatch+patch_index;
                else
                    error("Unknown patch type: " + loc(3:end) + ".")
                end
            end
        end


        %%% initialize block with bead indexed by ib_conn connected by r12_conn to r_source 
        function [b,failed,r_other] = init_positions(b,p,r_source,r12_uv_conn,r12_eq_conn,ib_conn,R,r_other)
            max_attempts = 100;

            %%% determine if connection direction should be random
            randomize_conn_dir = false;
            if r12_uv_conn == false
                randomize_conn_dir = true;
            end

            %%% determine if block direction should be random
            randomize_block_dir = false;
            if R == false
                randomize_block_dir = true;
            end

            %%% placement attempt loop
            attempts = 0;
            while true

                %%% check block attempts
                if attempts == max_attempts
                    failed = true;
                    return
                end

                %%% get position of connected bead
                if randomize_conn_dir
                    r12_uv_conn = ars.randUnitVec();
                end

                %%% get direction of block
                if randomize_block_dir
                    z_basis = ars.randUnitVec();
                    y_basis = ars.unitVector(cross(z_basis,ars.randUnitVec()));
                    x_basis = cross(y_basis,z_basis);
                    R = [x_basis,y_basis,z_basis];
                end
    
                %%% real bead absolute positions
                r_start = r_source + r12_eq_conn.*r12_uv_conn;
                com = r_start - R*b.r_internal(:,ib_conn);
                b.r = zeros(3,b.n);
                for ib = 1:b.n_real
                    b.r(:,ib) = com + R*b.r_internal(:,ib);
                    overlap = ars.checkOverlap(b.r(:,ib),r_other,p.r12_cut_WCA,p.dbox);
                    if overlap
                        break
                    end
                end

                %%% reset if overlap
                if overlap
                    attempts = attempts + 1;
                    continue
                end
    
                %%% patch bead absolute positions
                for pi = 1:b.npatch
                    ib = b.n-b.npatch+pi;
                    b.r(:,ib) = com + R*b.r_internal(:,ib);
                end

                %%% block successfully initiated
                r_other = ars.myHorzcat(r_other,b.r(:,1:b.n-b.npatch));
                failed = false;
                return
            end
        end
    
    end
    methods (Static)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Static Functions %%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%% get patterns (all lengths are in r12_helix units)
        function [pattern,hiHollow,hiL,hiR,hiM] = setpat(pattern_label)

            if pattern_label == "1HB"
                %%% all helices
                pattern = [0;...
                           0];
                %%% helix identification
                hiHollow = [];
                hiL = 1;
                hiR = 1;
                hiM = 1;

            elseif pattern_label == "sq8HB"
                %%% all helices
                pattern = [-1  0  1 -1  1 -1  0  1;...
                           -1 -1 -1  0  0  1  1  1];
                %%% helix identification
                hiHollow = [];
                hiL = 1;
                hiR = 3;
                hiM = 2;

            elseif pattern_label == "hex6HB"
                q = sqrt(3)/2;
                %%% all helices
                pattern = [-1.0 -0.5 0.5 1.0 0.5 -0.5;...
                           0    -q   -q  0   q   q  ];
                %%% helix identification
                hiHollow =  [];
                hiL = 1;
                hiR = 3;
                hiM = 2;

            elseif pattern_label == "hex16HB"
                q = sqrt(3)/2;
                %%% all helices
                pattern = [0 2 3 5  0.5 1.5 3.5 4.5  0.0 2.0 3.0 5.0  0.5 1.5 3.5 4.5;...
                           0 0 0 0  1*q 1*q 1*q 1*q  2*q 2*q 2*q 2*q  3*q 3*q 3*q 3*q];
                %%% helix identification
                hiHollow = [6,7,10,11];
                hiL = 1;
                hiR = 4;
                hiM = 2;

            else
                error("Unknown block pattern.")
            end
        end

    end
end