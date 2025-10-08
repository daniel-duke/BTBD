%%% origami class for BTBD
classdef origami
    properties
        rseed               % random seed for initialization
        rigid               % origami or block
        bs                  % block objects
        nconn               % number of connections
        conns_bis           % connection block indices
        conns_ibs           % connection bead indices
        conns_pot           % connection potential
        nangle              % number of angles
        angles_bis          % angle block indices
        angles_ibs          % angle bead indices
        angles_apot         % angle potential
        angles_theta_init   % angle initial theta
        nlink5              % number of 5p linkers
        link5s_index        % 5p linkers linker name
        link5s_io           % 5p linkers bead index within origami
        nlink3              % number of 3p linker ends
        link3s_index        % 3p linkers index
        link3s_io           % 3p linkers bead index within origami
        n                   % total number of beads
        r                   % total positions
        get_bi              % map io to bi
        get_ib              % map io to ib
        get_io              % map bi and ib to io
        index               % origami type
    end

    methods
        %%% constructor
        function o = origami(index)
            o.rseed = 0;
            o.rigid = "block";
            o.bs = [];
            o.nconn = 0;
            o.conns_bis = [];
            o.conns_ibs = [];
            o.conns_pot = bond_pot.empty;
            o.nangle = 0;
            o.angles_bis = [];
            o.angles_ibs = [];
            o.angles_apot = angle_pot.empty;
            o.angles_theta_init = [];
            o.nlink5 = 0;
            o.link5s_index = [];
            o.link5s_io = [];
            o.nlink3 = 0;
            o.link3s_index = [];
            o.link3s_io = [];
            o.index = index;
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
        function o = add_conn(o,bi1,loc1,bi2,loc2,pot)
            o.nconn = o.nconn + 1;
            o.conns_bis(:,o.nconn) = [bi1;bi2];
            ib1 = o.bs(bi1).interpret_loc(loc1);
            ib2 = o.bs(bi2).interpret_loc(loc2);
            o.conns_ibs(:,o.nconn) = [ib1;ib2];
            o.conns_pot(o.nconn) = pot;
        end


        %%% add angle to origami
        function o = add_angle(o,bi1,loc11,loc12,bi2,loc21,loc22,apot,theta_init)
            o.nangle = o.nangle + 1;
            o.angles_bis(:,o.nangle) = [bi1;bi1;bi2;bi2];
            ib11 = o.bs(bi1).interpret_loc(loc11);
            ib12 = o.bs(bi1).interpret_loc(loc12);
            ib21 = o.bs(bi2).interpret_loc(loc21);
            ib22 = o.bs(bi2).interpret_loc(loc22);
            o.angles_ibs(:,o.nangle) = [ib11;ib12;ib21;ib22];
            o.angles_apot(o.nangle) = apot;
            o.angles_theta_init(o.nangle) = theta_init;
        end


        %%% add linker to origami
        function o = add_linker(o,is_5p,li,bi,loc)
            if bi == "A"
                bi_min = 1;
                bi_max = length(o.bs);
            elseif bi == "B"
                bi_min = 1;
                bi_max = length(o.bs) - 1;
            else
                bi_min = str2double(bi);
                bi_max = str2double(bi);
            end
            for bi = bi_min:bi_max
                ib = o.bs(bi).interpret_loc(loc);
                io = o.get_io(bi,ib);
                if is_5p
                    o.nlink5 = o.nlink5 + 1;
                    o.link5s_index(o.nlink5) = li;
                    o.link5s_io(o.nlink5) = io;
                else
                    o.nlink3 = o.nlink3 + 1;
                    o.link3s_index(o.nlink3) = li;
                    o.link3s_io(o.nlink3) = io;
                end
            end
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


        %%% check if bead is patch
        function result = is_patch(o,io)
            b = o.bs(o.get_bi(io));
            result = o.get_ib(io) > b.n_real;
        end


        %%% ensure connections are not too strained
        function failed = are_conns_overstretched(o,p,bi_max)
            for ci = 1:o.nconn
                if o.conns_bis(1,ci) <= bi_max && o.conns_bis(2,ci) <= bi_max
                    r1 = o.bs(o.conns_bis(1,ci)).r(:,o.conns_ibs(1,ci));
                    r2 = o.bs(o.conns_bis(2,ci)).r(:,o.conns_ibs(2,ci));
                    U = o.conns_pot(ci).calc_energy(norm(r1-r2));
                    if U > p.U_overstretched
                        failed = true;
                        return
                    end
                end
            end
            failed = false;
        end


        %%% initialize origami positions
        function [o,failed,r_other] = init_positions(o,p,r_other)
            max_attempts_conf = 1000;
            max_attempts_rot = 1000;
            max_attempts_place = 1000;
            tolerance_rot = 0.001;

            %%% scale internal positions of structural beads to real units
            for bi = 1:length(o.bs)
                o.bs(bi).r_internal(:,1:o.bs(bi).n_real) = o.bs(bi).r_internal(:,1:o.bs(bi).n_real).*[p.r12_helix;p.r12_helix;p.r12_bead];
            end

            %%% configuraiton attempt loop
            disp("Looking for configuration...")
            attempts_conf = 0;
            while true

                %%% check configuration attempts
                if attempts_conf == max_attempts_conf
                    failed = true;
                    return
                end

                %%% set random seed
                if o.rseed ~= 0
                    attempts_conf = max_attempts_conf-1;
                    rstate = rng;
                    rng(o.rseed)
                end

                %%% initialize avoided positions
                r_other_origami = [];

                %%% add first block
                r_source = zeros(3,1); r12_eq_conn = 0; ib_conn = 1;
                [o.bs(1),~,r_other_origami] = o.bs(1).init_positions(p,r_source,false,r12_eq_conn,ib_conn,false,r_other_origami);

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
                                break
                            end
                        elseif o.conns_bis(1,ci) == bi
                            if o.conns_bis(2,ci) < bi
                                b0 = o.conns_bis(2,ci);
                                ib_b0_conn1 = o.conns_ibs(2,ci);
                                ib_b1_conn1 = o.conns_ibs(1,ci);
                                r12_eq_conn1 = o.conns_pot(ci).r12_eq;
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
                                ib_end_b1_angle1 = o.angles_ibs(4,ai);
                                break
                            end
                        end
                        if o.angles_bis(3,ai) == b0 && o.angles_ibs(3,ai) == ib_b0_conn1
                            if o.angles_bis(2,ai) == bi && o.angles_ibs(2,ai) == ib_b1_conn1
                                r12_uv_patches_b0_angle1 = ars.unitVector( o.bs(b0).r(:,o.angles_ibs(4,ai)) - o.bs(b0).r(:,o.angles_ibs(3,ai)) );
                                theta_init_angle1 = o.angles_theta_init(ai);
                                ib_end_b1_angle1 = o.angles_ibs(1,ai);
                                break
                            end
                        end
                    end

                    %%% if angle found, calculate block rotation
                    R = false;
                    if theta_init_angle1 ~= false

                        %%% look for second connection between the two blocks
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

                        %%% no second connection
                        if r12_eq_conn2 == false

                            %%% partially random conn1 direction (constrained by angle1)
                            r12_uv_perp = ars.unitVector(cross(r12_uv_patches_b0_angle1,ars.randUnitVec()));
                            r12_uv_conn1 = cosd(180-theta_init_angle1)*r12_uv_patches_b0_angle1 + sind(180-theta_init_angle1)*r12_uv_perp;

                            %%% partially random block direction (constrained by angle1)
                            r12_uv_perp = ars.unitVector(cross(r12_uv_conn1,ars.randUnitVec()));
                            r12_uv_patches_b1_angle1 = cosd(180-theta_init_angle1)*r12_uv_conn1 + sind(180-theta_init_angle1)*r12_uv_perp;
                            r12_uv_internal_patches_b1 = ars.unitVector(o.bs(bi).r_internal(:,ib_end_b1_angle1) - o.bs(bi).r_internal(:,ib_b1_conn1));
                            R = o.getRotation(r12_uv_patches_b1_angle1,r12_uv_internal_patches_b1,0);
                        
                        %%% match the second connection
                        else

                            %%% look for angle over second connection
                            theta_init_angle2 = false;
                            for ai = 1:o.nangle
                                if o.angles_bis(2,ai) == b0 && o.angles_ibs(2,ai) == ib_b0_conn2
                                    if o.angles_bis(3,ai) == bi && o.angles_ibs(3,ai) == ib_b1_conn2
                                        r12_uv_patches_b0_angle2 = ars.unitVector( o.bs(b0).r(:,o.angles_ibs(2,ai)) - o.bs(b0).r(:,o.angles_ibs(1,ai)) );
                                        theta_init_angle2 = o.angles_theta_init(ai);
                                        break
                                    end
                                end
                                if o.angles_bis(3,ai) == b0 && o.angles_ibs(3,ai) == ib_b0_conn1
                                    if o.angles_bis(2,ai) == bi && o.angles_ibs(2,ai) == ib_b1_conn1
                                        r12_uv_patches_b0_angle2 = ars.unitVector( o.bs(b0).r(:,o.angles_ibs(4,ai)) - o.bs(b0).r(:,o.angles_ibs(3,ai)) );
                                        theta_init_angle2 = o.angles_theta_init(ai);
                                        break
                                    end
                                end
                            end

                            %%% block rotation attempt loop
                            attempts_rot = 0;
                            while true

                                %%% check configuration attempts
                                if attempts_rot == max_attempts_rot
                                    failed = true;
                                    return
                                end

                                %%% connection direction
                                r12_uv_conn12_b0 = ars.unitVector(o.bs(b0).r(:,ib_b0_conn2) - o.bs(b0).r(:,ib_b0_conn1));
                                r12_uv_perp = (2*randi(2)-3)*ars.unitVector(cross(r12_uv_patches_b0_angle1,r12_uv_conn12_b0));
                                r12_uv_conn1 = cosd(180-theta_init_angle1)*r12_uv_patches_b0_angle1 + sind(180-theta_init_angle1)*r12_uv_perp;

                                %%% block direction
                                r12_uv_perp = (2*randi(2)-3)*ars.unitVector(cross(r12_uv_conn1,r12_uv_conn12_b0));
                                r12_uv_patches_b1_angle1 = cosd(180-theta_init_angle1)*r12_uv_conn1 + sind(180-theta_init_angle1)*r12_uv_perp;
                                r12_uv_internal_patches_b1 = ars.unitVector(o.bs(bi).r_internal(:,ib_end_b1_angle1) - o.bs(bi).r_internal(:,ib_b1_conn1));
                                [R,failed_rot,r12_conn2] = o.optimizeR(r12_uv_patches_b1_angle1,r12_uv_internal_patches_b1,b0,bi,ib_b0_conn1,ib_b1_conn1,ib_b0_conn2,ib_b1_conn2,r12_uv_conn1,r12_eq_conn1,r12_eq_conn2,tolerance_rot);
                                if failed_rot
                                    attempts_rot = attempts_rot + 1;
                                    continue
                                end

                                %%% check second angle
                                if theta_init_angle2 ~= 0
                                    if abs(dot(r12_uv_patches_b0_angle2, ars.unitVector(r12_conn2)) - cosd(180-theta_init_angle2)) > tolerance_rot
                                        attempts_rot = attempts_rot + 1;
                                        continue
                                    end
                                end

                                %%% block rotation found
                                break
                            end
                        end
                    end

                    %%% add block
                    r_source = o.bs(b0).r(:,ib_b0_conn1);
                    [o.bs(bi),failed_block,r_other_origami] = o.bs(bi).init_positions(p,r_source,r12_uv_conn1,r12_eq_conn1,ib_b1_conn1,R,r_other_origami);
                    if failed_block
                        break
                    end
                end

                %%% reset if internal block overlap
                if failed_block
                    attempts_conf = attempts_conf + 1;
                    continue
                end

                %%% reset if connection overstretch
                if are_conns_overstretched(o,p,bi)
                    attempts_conf = attempts_conf + 1;
                    continue
                end

                %%% origami configuration found
                break
            end

            %%% reset random seed
            if o.rseed ~= 0
                rng(rstate)
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
                    T = eye(3);
                    com = zeros(3,1);
                else
                    z_basis = ars.randUnitVec();
                    y_basis = ars.unitVector(cross(z_basis,ars.randUnitVec()));
                    x_basis = cross(y_basis,z_basis);
                    T = [x_basis,y_basis,z_basis];
                    com = (rand(3,1)-0.5)*p.dbox;
                end

                %%% rotate and move origami
                for io = 1:o.n
                    bi = o.get_bi(io);
                    ib = o.get_ib(io);
                    r_propose = com + T*o.bs(bi).r(:,ib);
                    propose(bi).r(:,ib) = r_propose;
                    if ib <= o.bs(bi).n_real
                        overlap = ars.checkOverlap(r_propose,r_other,p.r12_cut_WCA,p.dbox);
                        if overlap
                            break
                        end
                    end
                end

                %%% reset if overlap
                if overlap
                    attempts_place = attempts_place + 1;
                    continue
                end

                %%% origami successfully initiated
                o.r = [propose.r];
                for bi = 1:length(o.bs)
                    o.bs(bi).r = propose(bi).r;
                    r_real = o.bs(bi).r(:,1:o.bs(bi).n_real);
                    r_other = ars.myHorzcat(r_other,r_real);
                end
                failed = false;
                return
            end
        end


        %%% optimize rotation matrix around a connected angle to satisfy a second connection
        function [R,failed,r12_conn2] = optimizeR(o,r12_axis_b1,r12_axis_internal_b1,b0,b1,ib_b0_conn,ib_b1_conn,ib_b0_conn2,ib_b1_conn2,r12_uv_conn,r12_eq_conn,r12_eq_conn2,tolerance)

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
            finding_root = true;
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
                R = o.getRotation(r12_axis_b1,r12_axis_internal_b1,phi);
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

        function R = getRotation(v_f1,v_f2,phi)
            a = v_f1/norm(v_f1); 
            b = v_f2/norm(v_f2);
            ref = [0;0;1];
            if norm(ref-a) == 0 || norm(ref-b) == 0
                ref = [1;0;0];
                if norm(ref-a) == 0 || norm(ref-b)
                    ref = [0;1;0];
                end
            end
            ua = ref - a*(a'*ref); ua = ua/norm(ua); va = cross(a,ua);
            ub = ref - b*(b'*ref); ub = ub/norm(ub); vb = cross(b,ub);
            A = [a ua va]; 
            B = [b ub vb];
            R_0 = B*A';
            u = b; U = [0 -u(3) u(2); u(3) 0 -u(1); -u(2) u(1) 0];
            R_rot = eye(3) + sind(phi)*U + (1-cosd(phi))*(U*U);
            R = R_rot * R_0;
            R = R';
        end

    end
end