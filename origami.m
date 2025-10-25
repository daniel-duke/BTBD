%%% origami class for BTBD
classdef origami
    properties
        label               % origami name
        rigid               % origami or block
        bs                  % block objects
        nconn               % number of connections
        conns_bis           % connection block indices
        conns_ibs           % connection bead indices
        conns_pot           % connection potential
        conns_status        % connection status
        nangle              % number of angles
        angles_bis          % angle block indices
        angles_ibs          % angle bead indices
        angles_apot         % angle potential
        angles_theta_init   % angle initial theta
        angles_axis_init    % angle initial axis
        angles_status       % angle status
        ndihedral           % number of dihedrals
        dihedrals_bis       % dihedral block indices
        dihedrals_ibs       % dihedral bead indices
        dihedrals_dpot      % dihedral potential
        n                   % total number of beads
        r                   % total positions
        get_bi              % map io to bi
        get_ib              % map io to ib
        get_io              % map bi and ib to io
    end

    methods
        %%% constructor
        function o = origami(label)
            if nargin > 0
                o.label = label;
                o.rigid = "block";
                o.bs = [];
                o.nconn = 0;
                o.conns_bis = [];
                o.conns_ibs = [];
                o.conns_pot = bond_pot.empty;
                o.conns_status = [];
                o.nangle = 0;
                o.angles_bis = [];
                o.angles_ibs = [];
                o.angles_apot = angle_pot.empty;
                o.angles_theta_init = [];
                o.angles_axis_init = [];
                o.angles_status = [];
                o.ndihedral = 0;
                o.dihedrals_bis = [];
                o.dihedrals_ibs = [];
                o.dihedrals_dpot = dihedral_pot.empty;
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Utility Functions %%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%% add block to origami
        function o = add_block(o,b)
            o.bs = [o.bs b];
            o.n = sum([o.bs.n]);
            o.r = zeros(3,o.n);
            [o.get_bi,o.get_ib,o.get_io] = map_indices(o);
        end


        %%% add connection to origami
        function o = add_conn(o,bi1,patch1,bi2,patch2,pot)
            o.nconn = o.nconn + 1;
            o.conns_bis(:,o.nconn) = [bi1;bi2];
            ib1 = o.bs(bi1).get_ib_patch(patch1);
            ib2 = o.bs(bi2).get_ib_patch(patch2);
            o.conns_ibs(:,o.nconn) = [ib1;ib2];
            o.conns_pot(o.nconn) = pot;
            o.conns_status(o.nconn) = 1;
        end


        %%% add angle to origami
        function o = add_angle(o,bi1,patch11,patch12,bi2,patch21,patch22,apot)
            o.nangle = o.nangle + 1;
            o.angles_bis(:,o.nangle) = [bi1;bi1;bi2;bi2];
            ib11 = o.bs(bi1).get_ib_patch(patch11);
            ib12 = o.bs(bi1).get_ib_patch(patch12);
            ib21 = o.bs(bi2).get_ib_patch(patch21);
            ib22 = o.bs(bi2).get_ib_patch(patch22);
            o.angles_ibs(:,o.nangle) = [ib11;ib12;ib21;ib22];
            o.angles_apot(o.nangle) = apot;
            o.angles_theta_init(o.nangle) = apot.theta_eq;
            o.angles_axis_init(:,o.nangle) = [0;0;0];
            o.angles_status(o.nangle) = 1;
        end


        %%% add dihedral to origami
        function o = add_dihedral(o,bi1,patch1,bi2,patch2,bi3,patch3,bi4,patch4,dpot)
            o.ndihedral = o.ndihedral + 1;
            o.dihedrals_bis(:,o.ndihedral) = [bi1;bi2;bi3;bi4];
            ib1 = o.bs(bi1).get_ib_patch(patch1);
            ib2 = o.bs(bi2).get_ib_patch(patch2);
            ib3 = o.bs(bi3).get_ib_patch(patch3);
            ib4 = o.bs(bi4).get_ib_patch(patch4);
            o.dihedrals_ibs(:,o.ndihedral) = [ib1;ib2;ib3;ib4];
            o.dihedrals_dpot(o.ndihedral) = dpot;
        end


        %%% map between origami index and block indices
        function [get_bi,get_ib,get_io] = map_indices(o)
            get_bi = zeros(1,o.n);
            get_ib = zeros(1,o.n);
            get_io = zeros(length([o.bs]),max([o.bs.n]));
            io = 0;
            for bi = 1:length([o.bs])
                for ib = 1:o.bs(bi).n
                    io = io+1;
                    get_bi(io) = bi;
                    get_ib(io) = ib;
                    get_io(bi,ib) = io;
                end
            end
        end


        %%% ensure connections are not too strained
        function failed = are_conns_strained(o,p,bi)
            failed = false;
            for ci = 1:o.nconn
                if o.conns_bis(1,ci) <= bi && o.conns_bis(2,ci) <= bi
                    if o.conns_bis(1,ci) == bi || o.conns_bis(2,ci) == bi
                        r1 = o.bs(o.conns_bis(1,ci)).r(:,o.conns_ibs(1,ci));
                        r2 = o.bs(o.conns_bis(2,ci)).r(:,o.conns_ibs(2,ci));
                        U = o.conns_pot(ci).calc_energy(norm(r1-r2));
                        if U > p.U_strained
                            failed = true;
                            return
                        end
                    end
                end
            end
        end

        
        %%% ensure angles are not too strained
        function failed = are_angles_strained(o,p,bi)
            failed = false;
            for ai = 1:o.nangle
                if o.angles_bis(1,ai) <= bi && o.angles_bis(3,ai) <= bi
                    if o.angles_bis(1,ai) == bi || o.angles_bis(3,ai) == bi
                        r1 = o.bs(o.angles_bis(1,ai)).r(:,o.angles_ibs(1,ai));
                        r2 = o.bs(o.angles_bis(2,ai)).r(:,o.angles_ibs(2,ai));
                        r3 = o.bs(o.angles_bis(3,ai)).r(:,o.angles_ibs(3,ai));
                        r4 = o.bs(o.angles_bis(4,ai)).r(:,o.angles_ibs(4,ai));
                        r12_uv = ars.unitVector(r2-r1);
                        r23_uv = ars.unitVector(r3-r2);
                        r34_uv = ars.unitVector(r4-r3);
                        theta1 = 180-acosd(dot(r12_uv,r23_uv));
                        theta2 = 180-acosd(dot(r23_uv,r34_uv));
                        U1 = o.angles_apot(ai).calc_energy(theta1,o.angles_theta_init(ai));
                        U2 = o.angles_apot(ai).calc_energy(theta2,o.angles_theta_init(ai));
                        if U1 > p.U_strained || U2 > p.U_strained
                            failed = true;
                            return
                        end
                    end
                end
            end
        end


        %%% ensure dihedrals are not too strained
        function failed = are_dihedrals_strained(o,p,bi)
            failed = false;
            for di = 1:o.ndihedral
                if o.dihedrals_bis(1,di) <= bi && o.dihedrals_bis(2,di) <= bi && o.dihedrals_bis(3,di) <= bi && o.dihedrals_bis(4,di) <= bi
                    if o.dihedrals_bis(1,di) == bi || o.dihedrals_bis(2,di) == bi || o.dihedrals_bis(3,di) == bi || o.dihedrals_bis(4,di) == bi
                        r1 = o.bs(o.dihedrals_bis(1,di)).r(:,o.dihedrals_ibs(1,di));
                        r2 = o.bs(o.dihedrals_bis(2,di)).r(:,o.dihedrals_ibs(2,di));
                        r3 = o.bs(o.dihedrals_bis(3,di)).r(:,o.dihedrals_ibs(3,di));
                        r4 = o.bs(o.dihedrals_bis(4,di)).r(:,o.dihedrals_ibs(4,di));
                        r12_uv = ars.unitVector(r2-r1);
                        r23_uv = ars.unitVector(r3-r2);
                        r34_uv = ars.unitVector(r4-r3);
                        phi = ars.calcDihedral(r12_uv,r23_uv,r34_uv);
                        U = o.dihedrals_dpot(di).calc_energy(phi);
                        if U > p.U_strained
                            failed = true;
                            return
                        end
                    end
                end
            end
        end


        %%% initialize origami positions
        function [o,failed,r_other] = init_positions(o,p,r_other)
            max_attempts_block = 100;
            max_attempts_conf = 10;
            max_attempts_place = 100;
            tolerance = 0.001;

            %%% configuration attempt loop
            disp("Looking for configuration...")
            attempts_conf = 0;
            while true

                %%% check configuration attempts
                if attempts_conf == max_attempts_conf
                    failed = true;
                    return
                end

                %%% initialize avoided positions
                r_other_origami = [];

                %%% add first block
                ib_conn = 1; r_conn = zeros(3,1); R = ars.randBasis();
                [o.bs(1),~] = o.bs(1).init_positions(p,ib_conn,r_conn,R,r_other_origami);
                r_other_origami = o.bs(1).append_positions(r_other_origami);

                %%% loop over remaining blocks
                failed_block = false;
                for bi = 2:length(o.bs)
    
                    %%% find first connection between block and any previous block
                    for ci = 1:o.nconn
                        if o.conns_bis(2,ci) == bi
                            if o.conns_bis(1,ci) < bi
                                b0 = o.conns_bis(1,ci);
                                ib_b0_conn1 = o.conns_ibs(1,ci);
                                ib_b1_conn1 = o.conns_ibs(2,ci);
                                r12_eq_conn1 = o.conns_pot(ci).r12_eq;
                                R_b0_conn1 = o.bs(o.conns_bis(1,ci)).R;
                                break
                            end
                        elseif o.conns_bis(1,ci) == bi
                            if o.conns_bis(2,ci) < bi
                                b0 = o.conns_bis(2,ci);
                                ib_b0_conn1 = o.conns_ibs(2,ci);
                                ib_b1_conn1 = o.conns_ibs(1,ci);
                                r12_eq_conn1 = o.conns_pot(ci).r12_eq;
                                R_b0_conn1 = o.bs(o.conns_bis(2,ci)).R;
                                break
                            end
                        elseif ci == o.nconn
                            error("Unconnected block found")
                        end
                    end
    
                    %%% look for angle over first connection
                    theta_init_angle1 = false;
                    for ai = 1:o.nangle
                        if o.angles_bis(2,ai) == b0 && o.angles_ibs(2,ai) == ib_b0_conn1
                            if o.angles_bis(3,ai) == bi && o.angles_ibs(3,ai) == ib_b1_conn1
                                r12_uv_patches_b0_angle1 = ars.unitVector( o.bs(b0).r(:,o.angles_ibs(2,ai)) - o.bs(b0).r(:,o.angles_ibs(1,ai)) );
                                theta_init_angle1 = o.angles_theta_init(ai);
                                axis_init_angle1 = o.angles_axis_init(:,ai);
                                ib_end_b1_angle1 = o.angles_ibs(4,ai);
                                break
                            end
                        end
                        if o.angles_bis(3,ai) == b0 && o.angles_ibs(3,ai) == ib_b0_conn1
                            if o.angles_bis(2,ai) == bi && o.angles_ibs(2,ai) == ib_b1_conn1
                                r12_uv_patches_b0_angle1 = ars.unitVector( o.bs(b0).r(:,o.angles_ibs(4,ai)) - o.bs(b0).r(:,o.angles_ibs(3,ai)) );
                                theta_init_angle1 = o.angles_theta_init(ai);
                                axis_init_angle1 = o.angles_axis_init(:,ai);
                                ib_end_b1_angle1 = o.angles_ibs(1,ai);
                                break
                            end
                        end
                    end

                    %%% look for second connection between the two blocks
                    if theta_init_angle1 ~= 0
                        r12_eq_conn2 = false;
                        for ci = ci+1:o.nconn
                            if o.conns_bis(2,ci) == bi
                                if o.conns_bis(1,ci) == b0
                                    ib_b0_conn2 = o.conns_ibs(1,ci);
                                    ib_b1_conn2 = o.conns_ibs(2,ci);
                                    r12_eq_conn2 = o.conns_pot(ci).r12_eq;
                                    break
                                end
                            elseif o.conns_bis(1,ci) == bi
                                if o.conns_bis(2,ci) == b0
                                    ib_b0_conn2 = o.conns_ibs(2,ci);
                                    ib_b1_conn2 = o.conns_ibs(1,ci);
                                    r12_eq_conn2 = o.conns_pot(ci).r12_eq;
                                    break
                                end
                            end
                        end
                    end

                    %%% block attempt loop
                    attempts_block = 0;
                    while true

                        %%% check block attempts
                        if attempts_block == max_attempts_block
                            failed_block = true;
                            break
                        end
                        
                        %%% no angle over connection
                        if theta_init_angle1 == 0

                            %%% connection position
                            r_b1_conn1 = o.bs(b0).r(:,ib_b0_conn1) + r12_eq_conn1*ars.randUnitVec();
    
                            %%% block rotation
                            R = ars.randBasis();

                        %%% angle over connection
                        else

                            %%% bending axis
                            if all(axis_init_angle1==0)
                                r12_axis = cross(r12_uv_patches_b0_angle1,ars.randUnitVec());
                            else
                                r12_axis = R_b0_conn1*axis_init_angle1;
                                if abs(dot(r12_uv_patches_b0_angle1,r12_axis)) > tolerance
                                    error("Angle axis not perpendicular to first bond in angle.")
                                end
                            end

                            %%% connection position
                            r12_uv_bend_dir = ars.unitVector(cross(r12_uv_patches_b0_angle1,r12_axis));
                            r12_uv_conn1 = cosd(180-theta_init_angle1)*r12_uv_patches_b0_angle1 + sind(180-theta_init_angle1)*r12_uv_bend_dir;
                            r_b1_conn1 = o.bs(b0).r(:,ib_b0_conn1) + r12_eq_conn1*r12_uv_conn1;

                            %%% block direction
                            r12_uv_bend_dir = ars.unitVector(cross(r12_uv_conn1,r12_axis));
                            r12_uv_patches_b1_angle1 = cosd(180-theta_init_angle1)*r12_uv_conn1 + sind(180-theta_init_angle1)*r12_uv_bend_dir;
                            r12_uv_internal_patches_b1 = ars.unitVector(o.bs(bi).r_internal(:,ib_end_b1_angle1) - o.bs(bi).r_internal(:,ib_b1_conn1));

                            %%% block rotation
                            if r12_eq_conn2 == 0
                                phi = rand*360;
                                R = o.getRotation(r12_uv_internal_patches_b1,r12_uv_patches_b1_angle1,phi);
                            else
                                [R,failed_block_attempt] = o.optimizeR(r12_uv_patches_b1_angle1,r12_uv_internal_patches_b1,b0,bi,ib_b0_conn1,ib_b1_conn1,ib_b0_conn2,ib_b1_conn2,r12_uv_conn1,r12_eq_conn1,r12_eq_conn2,tolerance);
                                if failed_block_attempt
                                    attempts_block = attempts_block + 1;
                                    continue
                                end
                            end
                        end

                        %%% place block
                        [o.bs(bi),failed_block_attempt] = o.bs(bi).init_positions(p,ib_b1_conn1,r_b1_conn1,R,r_other_origami);
                        if failed_block_attempt
                            attempts_block = attempts_block + 1;
                            continue
                        end

                        %%% check connections and angles and dihedrals
                        if are_conns_strained(o,p,bi) || are_angles_strained(o,p,bi) || are_dihedrals_strained(o,p,bi)
                            attempts_block = attempts_block + 1;
                            continue
                        end

                        %%% block successfully placed
                        r_other_origami = o.bs(bi).append_positions(r_other_origami);
                        failed_block = false;
                        break
                    end

                    %%% stop looping through blocks if one failed
                    if failed_block == true
                        break
                    end
                end

                %%% check if block loop succeeded
                if failed_block == true
                    attempts_conf = attempts_conf + 1;
                    continue
                end

                %%% origami configuration found
                break
            end

            %%% center positions
            r_real = [];
            for bi = 1:length(o.bs)
                r_real_block = o.bs(bi).r(:,1:o.bs(bi).n_real);
                r_real = ars.myHorzcat(r_real, r_real_block);
            end
            com_real = ars.calcCOM(r_real, p.dbox);
            for bi = 1:length(o.bs)
                o.bs(bi).r = o.bs(bi).r - com_real;
            end

            %%% placement attempt loop
            disp("Attempting to place...")
            attempts_place = 0;
            while true

                %%% check placement attempts
                if attempts_place == max_attempts_place
                    failed = true;
                    return
                end

                %%% initialize proposed positions
                propose = repmat(struct(),length(o.bs),1);
                for bi = 1:length(o.bs)
                    propose(bi).r = o.bs(bi).r;
                end

                %%% get random basis and position
                if attempts_place == 0
                    R = eye(3);
                    com = zeros(3,1);
                else
                    R = ars.randBasis();
                    com = (rand(3,1)-0.5)*p.dbox;
                end

                %%% rotate and move origami
                for io = 1:o.n
                    bi = o.get_bi(io);
                    ib = o.get_ib(io);
                    r_propose = com + R*o.bs(bi).r(:,ib);
                    propose(bi).r(:,ib) = r_propose;
                    if ib <= o.bs(bi).n_real
                        overlap = ars.checkOverlap(r_propose,r_other,p.r12_cut_WCA,p.dbox);
                        if overlap
                            break
                        end
                    end
                end

                %%% check if placement succeeded
                if overlap
                    attempts_place = attempts_place + 1;
                    continue
                end

                %%% update origami data
                o.r = [propose.r];
                for bi = 1:length(o.bs)
                    o.bs(bi).r = propose(bi).r;
                    r_real = o.bs(bi).r(:,1:o.bs(bi).n_real);
                    r_other = ars.myHorzcat(r_other,r_real);
                end

                %%% origami successfully placed
                failed = false;
                return
            end
        end


        %%% optimize rotation matrix around a connected angle to satisfy a second connection
        function [R,failed] = optimizeR(o,r12_axis_b1,r12_axis_internal_b1,b0,b1,ib_b0_conn,ib_b1_conn,ib_b0_conn2,ib_b1_conn2,r12_uv_conn,r12_eq_conn,r12_eq_conn2,tolerance)

            %%% notes
            % optimization starts by assuming the second connection
            % distance is always greater than or equal to the equilibrium
            % distance, and accordingly it attempts to minimize the
            % separation; but if it encounters a separation less than the
            % equilibrium distance, it starts searching for the exact
            % solution (essentially a root finding problem).

            %%% parameters
            max_iterations = 1000;
            initializing = true;
            finding_root = false;
            phi = rand*360;
            phi_prev = phi;
            step = 10;

            %%% numerical solver loop
            iteration = 0;
            while true

                %%% should never happen
                if iteration == max_iterations
                    failed = true;
                    break
                end

                %%% calculate distance of second connection
                R = o.getRotation(r12_axis_internal_b1,r12_axis_b1,phi);
                r_conn_bi = o.bs(b0).r(:,ib_b0_conn) + r12_uv_conn.*r12_eq_conn;
                r_conn_bi_internal = o.bs(b1).r_internal(:,ib_b1_conn);
                r_conn2_bi_internal = o.bs(b1).r_internal(:,ib_b1_conn2);
                r_conn2_bi = r_conn_bi + R*(r_conn2_bi_internal-r_conn_bi_internal);
                r_conn2_b0 = o.bs(b0).r(:,ib_b0_conn2);
                r12_conn2 = r_conn2_bi-r_conn2_b0;
                d = norm(r12_conn2);

                %%% first two iterations
                if initializing == true
                    if phi - phi_prev == step
                        m = (d-d_prev)/step;
                        if m > 0
                            step = -step;
                        end
                        m_prev = m;
                        initializing = false;
                    end

                %%% searching for minimum distance
                elseif finding_root == false

                    %%% check for convergence
                    if abs(m) < tolerance
                        failed = false;
                        break

                    %%% check distance less than equilibrium
                    elseif d-r12_eq_conn2 < 0
                        finding_root = true;
                        d = d_prev;
                        phi = phi - step;
                    
                    %%% determine step
                    else
                        m = (d-d_prev)/step;
                        if sign(m) ~= sign(m_prev)
                            step = -step/2;
                        end
                        m_prev = m;
                    end

                %%% searching for exact distance
                else

                    %%% check for convergence
                    if abs(d-r12_eq_conn2) < tolerance
                        failed = false;
                        break

                    %%% determine step
                    elseif sign(d-r12_eq_conn2) ~= sign(d_prev-r12_eq_conn2)
                        step = -step/2;
                    end
                end

                %%% take step
                d_prev = d;
                phi = phi + step;
                iteration = iteration + 1;
            end
        end
        
    end
    methods (Static)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Static Functions %%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%% calculate the rotation matrix between frames
        function R = getRotation(v_f1,v_f2,phi)
            
            %%% notes
            % if the same vector is defined in two frames, the
            % transformation matrix calculated by this function rotates
            % the vector from frame 1 into frame 2; however, this does not
            % fully define the transformation, so an additional rotational
            % angle phi must be given.

            %%% normalize
            a = ars.unitVector(v_f1); 
            b = ars.unitVector(v_f2);

            %%% get reference
            ref = [0;0;1];
            if norm(ref-a) == 0 || norm(ref-b) == 0
                ref = [1;0;0];
                if norm(ref-a) == 0 || norm(ref-b)
                    ref = [0;1;0];
                end
            end

            %%% linear algebra
            ua = ars.unitVector(ref - a*(a'*ref));
            ub = ars.unitVector(ref - b*(b'*ref));
            va = cross(a,ua);
            vb = cross(b,ub);
            A = [a ua va]; 
            B = [b ub vb];
            R_0 = B*A';

            %%% rotation
            u = b; U = [0 -u(3) u(2); u(3) 0 -u(1); -u(2) u(1) 0];
            R_rot = eye(3) + sind(phi)*U + (1-cosd(phi))*(U*U);
            R = R_rot * R_0;
        end

    end
end