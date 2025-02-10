%%% tether class for BTBD
classdef tether
    properties
        n           %number of beads
        r           %positions
    end

    methods
        %%% constructor
        function t = tether(n)
            if nargin > 0
                t.n = n;
                t.r = zeros(3,n);
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%% initialize tether positions with first bead at r_start and next
        %%% bead pointing in direction dir
        function [r,overlap] = init(t,p,r_start,dir,r_other)
            r = zeros(3,t.n);

            %%% place first bead
            r(:,1) = r_start;
            overlap = ars.check_overlap(r(:,1),r_other,p.sigma,p.dbox);
            if overlap == true
                return
            end

            %%% place second bead
            r(:,2) = r(:,1) + p.r12_eq_tether.*ars.unit_vector(dir);
            overlap = ars.check_overlap(r(:,2),[r_other,r(:,1)],p.sigma,p.dbox);
            if overlap == true
                return
            end

            %%% place remaining beads
            max_attempts = 10;
            for attempts = 1:max_attempts
                for i = 3:t.n
                    r12 = ars.applyPBC(r(:,i-1) - r(:,i-2),p.dbox);
                    r(:,i) = ars.applyPBC(r(:,i-1) + addPartSH(r12,p.r12_eq_tether,p.k_x_tether,p.k_theta,p.kBT), p.dbox);
                    overlap = ars.ars.check_overlap(r(:,i),[r_other,r(:,1:i-1)],p.sigma,p.dbox);
                    if overlap == true
                        break
                    end
                end
                if overlap == false
                    return
                end
            end
        end

        %%% create random vector constrained by bending potential
        function r23 = addPartSH(r12,r12_eq,k_x,k_theta,kBT)
            phi = rand*2*pi;
            theta = sqrt(2*kBT/k_theta*log(1/rand));
            q = [cos(phi/2),sin(phi/2).*ars.unit_vector(r12)'];
            A = ars.quat_rot(q(1),q(2),q(3),q(4));
            perp_unit = ars.unit_vector(cross(r12,box_muller));
            r23_mag = r12_eq + (2*randi([0 1])-1)*sqrt(2*kBT/k_x*log(1/rand));
            r23 = r23_mag.*( cos(theta).*ars.unit_vector(r12) + sin(theta)*(A*perp_unit) );
        end

    end
end