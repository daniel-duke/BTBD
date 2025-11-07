%%% block class for BTBD
classdef block
    properties
        pattern         % helix locations in block xy-plane
        hi_hollow       % helix indices of hollow helices
        is_linear       % whether the pattern is 1HB
        nhelix          % number of helices
        n_helix         % number of beads in each helix
        r_internal      % internal positions (in bead units)
        R               % rotation matrix (block to lab coords)
        npatch          % number of patches
        patches         % patch name
        n_real          % number of real (structural) beads
        n               % number of beads (structural and patches)
        r               % bead positions
        omit            % whether to omit from simulation, for each bead
        hide            % whether to hide in visuals, for each bead
        get_hi          % map ib to hi
        get_zi          % map ib to zi
        get_ib          % map hi and zi to ib
    end

    methods
        %%% constructor
        function b = block(pattern_label,n_helix,p)
            if nargin > 0
                [b.pattern,b.hi_hollow,b.is_linear] = block.init(pattern_label);
                b.nhelix = size(b.pattern,2);
                b.n_helix = n_helix;
                b.R = eye(3);
                b.npatch = 0;
                b.patches = {};
                b.n_real = b.calc_nbead;
                b.n = b.n_real;
                b.r = zeros(3,b.n);
                b.omit = zeros(1,b.n);
                b.hide = zeros(1,b.n);
                [b.get_hi,b.get_zi,b.get_ib] = map_indices(b);
                b.r_internal = init_positions_internal(b,p);
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
            b.omit(b.n) = 0;
            b.hide(b.n) = 0;
            b.r_internal(:,b.n) = [x;y;z];
        end


        %%% calcualte number of structural beads in block
        function nbead = calc_nbead(b)
            nbead = size(b.pattern,2)*min(b.n_helix,2) + (size(b.pattern,2)-length(b.hi_hollow))*max(b.n_helix-2,0);
        end


        %%% get block index from geometrical indices
        function [get_hi,get_zi,get_ib] = map_indices(b)
            get_hi = zeros(1,b.n);
            get_zi = zeros(1,b.n);
            get_ib = zeros(b.nhelix,b.n_helix);
            for hi = 1:b.nhelix
                for zi = 1:b.n_helix
                    if zi ~= 1 && zi ~= b.n_helix && any(b.hi_hollow==hi)
                        ib = 0;
                    else
                        ib = size(b.pattern,2)*min(zi-1,1) + ...
                                (size(b.pattern,2)-length(b.hi_hollow))*max(zi-2,0) + hi;
                        if zi ~= 1 && zi ~= b.n_helix
                            ib = ib - sum(b.hi_hollow<hi);
                        end
                    end
                    get_hi(ib) = hi;
                    get_zi(ib) = zi;
                    get_ib(hi,zi) = ib;
                end
            end
        end


        %%% get block index from patch name
        function ib = get_ib_patch(b,patch)
            patch_index = find(strcmp(b.patches, patch));
            if ~isempty(patch_index)
                ib = b.n_real+patch_index;
            else
                error("Unknown patch type: " + patch + ".")
            end
        end


        %%% set the positions of the block within its own coordinate system
        function r_internal = init_positions_internal(b,p)
            r_internal = zeros(3,b.n);
            for ib = 1:b.n_real
                hi = b.get_hi(ib);
                zi = b.get_zi(ib);
                r_internal(1:2,ib) = p.r12_helix.*b.pattern(:,hi);
                r_internal(3,ib) = p.r12_bead.*(zi-1);
            end
        end

        
        %%% add real block positions to list of positions
        function r = append_positions(b,r)
            r = ars.myHorzcat(r,b.r(:,1:b.n-b.npatch));
        end


        %%% initialize block with bead indexed by ib_conn at r_conn
        function [b,failed] = init_positions(b,p,ib_conn,r_conn,R,r_other)

            %%% real bead absolute positions
            com = r_conn - R*b.r_internal(:,ib_conn);
            b.r = zeros(3,b.n);
            for ib = 1:b.n_real
                b.r(:,ib) = com + R*b.r_internal(:,ib);
                overlap = ars.checkOverlap(b.r(:,ib),r_other,p.r12_cut_WCA,p.dbox);
                if overlap
                    failed = true;
                    return
                end
            end

            %%% patch bead absolute positions
            for pi = 1:b.npatch
                ib = b.n-b.npatch+pi;
                b.r(:,ib) = com + R*b.r_internal(:,ib);
            end

            %%% block successfully initiated
            b.R = R;
            failed = false;
            return
        end
    
    end
    methods (Static)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Static Functions %%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%% get patterns
        function [pattern,hi_hollow,is_linear] = init(pattern_label)

            %%% no hollow helices by default
            hi_hollow = [];
            is_linear = 0;

            %%% block pattern
            switch pattern_label

                %%% single bundle
                case "1HB"
                    is_linear = 1;
                    pattern  = [ 0 ;...
                                 0 ];

                %%% square 4-helix bundle
                case "sq4HB"
                    h = 0.5;
                    pattern  = [ h  -h  -h   h ;...
                                 h   h  -h  -h ];

                %%% square 8-helix bundle
                case "sq8HB"
                    pattern  = [ 1  1  0  -1  -1  -1   0   1 ;...
                                 0  1  1   1   0  -1  -1  -1 ];

                %%% hexagonal 6-helix bundle
                case "hex6HB"
                    h = 0.5; q = sqrt(3)/2;
                    pattern  = [ 1  h  -h  -1  -h   h ;...
                                 0  q   q   0  -q  -q ];

                %%% error
                otherwise
                    error("Unknown block pattern.")
            end
        end

    end
end