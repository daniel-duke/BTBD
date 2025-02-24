%%% block class for BTBD
classdef block
    properties
        pattern     %helix locations in block xy-plane
        hiHollow    %helix indices of hollow helices
        hiL         %helix index for left connections
        hiR         %helix index for right connections
        hiM         %helix index for middle connections
        patches     %names of patches
        npatch      %number of patches
        r12_pat_pol %patch internal positions (theta,radius,z)
        r12_cart    %internal positions of all beads (x,y,z)
        n_xy        %number of helices
        n_z         %number of beads in each helix
        n_r         %number of real beads
        n           %number of beads in total
        r           %bead positions
    end

    methods
        %%% constructor
        function b = block(pattern_label,n_z)
            if nargin > 0
                [pattern,hiHollow,hiL,hiR,hiM] = block.setpat(pattern_label);
                b.pattern = pattern;
                b.hiHollow = hiHollow;
                b.hiL = hiL;
                b.hiR = hiR;
                b.hiM = hiM;
                b.patches = {};
                b.npatch = 0;
                b.r12_pat_pol = [];
                b.r12_cart = [];
                b.n_xy = size(b.pattern,2);
                b.n_z = n_z;
                b.n_r = size(b.pattern,2)*min(n_z,2) + (size(b.pattern,2)-length(b.hiHollow))*max(n_z-2,0);
                b.n = b.n_r;
                b.r = zeros(3,b.n);
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
        %%% add patch to block
        function b = add_patch(b,name,theta,radius,z)
            b.patches{end+1} = name;
            b.n = b.n + 1;
            b.r = [b.r zeros(3,1)];
            b.r12_pat_pol = [b.r12_pat_pol [theta;radius;z]];
            b.r12_cart = [b.r12_cart zeros(3,1)];
            b.npatch = b.npatch + 1;
        end


        %%% get block index from geometrical indices
        function index = get_i(b,hi,zi)
            if zi ~= 1 && zi ~= b.n_z && any(b.hiHollow==hi)
                index = 0;
            else
                index = size(b.pattern,2)*min(zi-1,1) + ...
                        (size(b.pattern,2)-length(b.hiHollow))*max(zi-2,0) + hi;
                if zi ~= 1 && zi ~= b.n_z
                    index = index - sum(b.hiHollow<hi);
                end
            end
        end


        %%% get block index from location
        function index = interpret_loc(b,loc)
            if loc(1) == 'P'
                patch_index = find(strcmp(b.patches, loc(3:end)));
                if ~isempty(patch_index)
                    index = b.n-b.npatch+patch_index;
                else
                    error("Unknown patch type: " + loc(3:end))
                end
            elseif loc(1) == 'B'
                if loc(3) == 'L'
                    hi = b.hiL;
                elseif loc(3) == 'M'
                    hi = b.hiM;
                elseif loc(3) == 'R'
                    hi = b.hiR;
                else
                    error("Unknown helix type: " + loc(3))
                end
                zi = str2double(loc(4:end));
                index = b.get_i(hi,zi);
            else
                error("Unknown connection type:" + loc(1))
            end
        end


        %%% initialize block with bead indexed by i_conn at r_conn and 
        %%% orientation defined by xy-plane normal vector (z-direction)
        function [b,overlap] = init(b,p,i_conn,r_conn,xy_normal,r_other)

            %%% set real bead internal positions
            for zi = 1:b.n_z
                for hi = 1:b.n_xy
                    index = b.get_i(hi,zi);
                    if index ~= 0
                        r12 = zeros(3,1);
                        r12(1:2) = p.r12_adj_block*(b.pattern(:,hi)-b.pattern(:,1));
                        r12(3) = p.r12_eq_block*(zi-1);
                        b.r12_cart(:,index) = r12;
                    end
                end
            end

            %%% set patch internal positions
            for i = 1:b.npatch
                index = b.n_r+i;
                theta = b.r12_pat_pol(1,i);
                radius = b.r12_pat_pol(2,i);
                z = b.r12_pat_pol(3,i);
                [x,y,z] = pol2cart(deg2rad(theta),radius*p.r12_adj_block,z*p.r12_eq_block);
                b.r12_cart(:,index) = [x;y;z];
            end

            %%% set real bead positions
            z_basis = ars.unitVector(xy_normal);
            y_basis = ars.unitVector(cross(z_basis,box_muller));
            x_basis = cross(y_basis,z_basis);
            T = [x_basis,y_basis,z_basis];
            origin = r_conn - T*b.r12_cart(:,i_conn);
            b.r = zeros(3,b.n);
            count = 0;
            for zi = 1:b.n_z
                for hi = 1:b.n_xy
                    index = b.get_i(hi,zi);
                    if index ~= 0
                        count = count+1;
                        b.r(:,index) = ars.applyPBC(origin + T*b.r12_cart(:,index), p.dbox);
                        overlap = ars.checkOverlap(b.r(:,index),r_other,p.sigma,p.dbox);
                        if overlap == true
                            return
                        end
                    end
                end
            end

            %%% set patch positions
            for i = 1:b.npatch
                index = b.n_r+i;
                b.r(:,index) = ars.applyPBC(origin + T*b.r12_cart(:,index), p.dbox);
            end
        end

    end
    methods (Static)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Static Functions %%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%% get patterns (all lengths are in r12_adj_block units)
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
                pattern = [0 1 2 0 2 0 1 2;...
                           0 0 0 1 1 2 2 2];
                %%% helix identification
                hiHollow = [];
                hiL = 1;
                hiR = 3;
                hiM = 2;

            elseif pattern_label == "hex6HB"
                q = sqrt(3)/2;
                %%% all helices
                pattern = [0.0 0.5 1.5 2.0 1.5 0.5;...
                           0   -q  -q  0   q   q  ];
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
                error("Unknown block pattern")
            end
        end

    end
end