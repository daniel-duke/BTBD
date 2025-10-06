%%% origami class for BTBD
classdef origami
    properties
        name                % origami name
        rigid               % origami or block
        bs                  % block objects
        nconn               % number of connections
        conns_bis           % connection block indices
        conns_ibs           % connection bead indices
        conns_pot           % name of connection potential
        conns_r12_eq        % connection equilibrium separation
        nangle              % number of angles
        angles_bis          % angle block indices
        angles_ibs          % angle bead indices
        angles_theta_eq     % angle equilibrium theta
        angles_theta_init   % angle initial theta
        nlink5              % number of 5p linkers
        link5s_name         % 5p linkers linker name
        link5s_io           % 5p linkers bead index within origami
        nlink3              % number of 3p linker ends
        link3s_name         % 3p linkers name
        link3s_io           % 3p linkers bead index within origami
        n                   % total number of beads
        r                   % total positions
        get_bi              % map io to bi
        get_ib              % map io to ib
        get_io              % map bi and ib to io
    end

    methods
        %%% constructor
        function o = origami(name)
            o.name = name;
            o.rigid = "block";
            o.bs = [];
            o.nconn = 0;
            o.conns_bis = [];
            o.conns_ibs = [];
            o.conns_pot = strings();
            o.conns_r12_eq = [];
            o.nangle = 0;
            o.angles_bis = [];
            o.angles_ibs = [];
            o.angles_theta_eq = [];
            o.angles_theta_init = [];
            o.nlink5 = 0;
            o.link5s_name = strings();
            o.link5s_io = [];
            o.nlink3 = 0;
            o.link3s_name = strings();
            o.link3s_io = [];
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
        function o = add_conn(o,bi1,loc1,bi2,loc2,pot,r12_eq)
            o.nconn = o.nconn + 1;
            o.conns_bis(:,o.nconn) = [bi1;bi2];
            ib1 = o.bs(bi1).interpret_loc(loc1);
            ib2 = o.bs(bi2).interpret_loc(loc2);
            o.conns_ibs(:,o.nconn) = [ib1;ib2];
            o.conns_pot(o.nconn) = pot;
            o.conns_r12_eq(o.nconn) = r12_eq;
        end


        %%% add angle to origami
        function o = add_angle(o,bi1,loc1,bi2,loc2,bi3,loc3,bi4,loc4,theta_eq,theta_init)
            o.nangle = o.nangle + 1;
            o.angles_bis(:,o.nangle) = [bi1;bi2;bi3;bi4];
            ib1 = o.bs(bi1).interpret_loc(loc1);
            ib2 = o.bs(bi2).interpret_loc(loc2);
            ib3 = o.bs(bi3).interpret_loc(loc3);
            ib4 = o.bs(bi4).interpret_loc(loc4);
            o.angles_ibs(:,o.nangle) = [ib1;ib2;ib3;ib4];
            o.angles_theta_eq(o.nangle) = theta_eq;
            o.angles_theta_init(o.nangle) = theta_init;
        end


        %%% add connection to origami
        function o = add_linker(o,name,is_5p,bi,loc)
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
                    o.link5s_name(o.nlink5) = name;
                    o.link5s_io(o.nlink5) = io;
                else
                    o.nlink3 = o.nlink3 + 1;
                    o.link3s_name(o.nlink3) = name;
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
        function failed = are_conns_overstretched(o,p,U_overstretched)
            for ci = 1:o.nconn
                r12_eq = o.conns_r12_eq(ci);
                r12_overstretched = r12_eq + sqrt(2*U_overstretched/p.k_x_conn);
                r1 = o.bs(o.conns_bis(1,ci)).r(:,o.conns_ibs(1,ci));
                r2 = o.bs(o.conns_bis(2,ci)).r(:,o.conns_ibs(2,ci));
                if norm(r1-r2) > r12_overstretched
                    failed = true;
                    return
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
            U_overstretched = 10;

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

                %%% initialize avoided positions
                r_other_origami = [];

                %%% add first block
                r_source = zeros(3,1);
                [o.bs(1),~,r_other_origami] = o.bs(1).init_positions(p,r_source,0,1,0,r_other_origami);

                %%% loop over remaining blocks
                failed_block = false;
                for bi = 2:length(o.bs)

                    %%% find first connection between block and any previous block
                    for ci = 1:o.nconn
                        if o.conns_bis(2,ci) == bi
                            if o.conns_bis(1,ci) < bi
                                b0 = o.conns_bis(1,ci);
                                ib_conn_b0 = o.conns_ibs(1,ci);
                                ib_conn_b1 = o.conns_ibs(2,ci);
                                r12_conn = o.conns_r12_eq(ci);
                                r12_conn_mag = r12_conn;
                                break
                            end
                        elseif o.conns_bis(1,ci) == bi
                            if o.conns_bis(2,ci) < bi
                                b0 = o.conns_bis(2,ci);
                                ib_conn_b0 = o.conns_ibs(2,ci);
                                ib_conn_b1 = o.conns_ibs(1,ci);
                                r12_conn = o.conns_r12_eq(ci);
                                r12_conn_mag = r12_conn;
                                break
                            end
                        elseif ci == o.nconn
                            error("Unconnected block found")
                        end
                    end

                    %%% look for angle over first connection
                    theta_init = 0;
                    for ai = 1:o.nangle
                        if o.angles_bis(2,ai) == b0 && o.angles_ibs(2,ai) == ib_conn_b0
                            if o.angles_bis(3,ai) == bi && o.angles_ibs(3,ai) == ib_conn_b1
                                r12_patches_b0 = o.bs(b0).r(:,o.angles_ibs(2,ai)) - o.bs(b0).r(:,o.angles_ibs(1,ai));
                                theta_init = o.angles_theta_init(ai);
                                ib_end_b1 = o.angles_ibs(4,ai);
                                break
                            end
                        end
                        if o.angles_bis(3,ai) == b0 && o.angles_ibs(3,ai) == ib_conn_b0
                            if o.angles_bis(2,ai) == bi && o.angles_ibs(2,ai) == ib_conn_b1
                                r12_patches_b0 = o.bs(b0).r(:,o.angles_ibs(4,ai)) - o.bs(b0).r(:,o.angles_ibs(3,ai));
                                theta_init = o.angles_theta_init(ai);
                                ib_end_b1 = o.angles_ibs(1,ai);
                                break
                            end
                        end
                    end

                    %%% if angle found, calculate directions
                    R = 0;
                    if theta_init > 0

                        %%% look for second connection between the two blocks
                        r12_eq_conn2 = 0;
                        for ci = ci+1:o.nconn
                            if o.conns_bis(2,ci) == bi
                                if o.conns_bis(1,ci) == b0
                                    ib_conn2_b0 = o.conns_ibs(1,ci);
                                    ib_conn2_b1 = o.conns_ibs(2,ci);
                                    r12_eq_conn2 = o.conns_r12_eq(ci);
                                    break
                                end
                            elseif o.conns_bis(1,ci) == bi
                                if o.conns_bis(2,ci) == b0
                                    ib_conn2_b0 = o.conns_ibs(2,ci);
                                    ib_conn2_b1 = o.conns_ibs(1,ci);
                                    r12_eq_conn2 = o.conns_r12_eq(ci);
                                    break
                                end
                            elseif ci == o.nconn
                                error("Unconnected block found")
                            end
                        end

                        %%% no second connection
                        if r12_eq_conn2 == 0

                            %%% connection direction
                            r12_perp = ars.unitVector(cross(r12_patches_b0,ars.randUnitVec()));
                            r12_conn = r12_conn_mag*(cosd(180-theta_init)*ars.unitVector(r12_patches_b0) + sind(180-theta_init)*r12_perp);

                            %%% block direction
                            r12_perp = ars.unitVector(cross(r12_conn,ars.randUnitVec()));
                            r12_patches_bi = cosd(180-theta_init)*ars.unitVector(r12_conn) + sind(180-theta_init)*r12_perp;
                            r12_patches_bi_internal = ars.unitVector(o.bs(bi).r_internal(:,ib_end_b1) - o.bs(bi).r_internal(:,ib_conn_b1));
                            R = o.getRotation(r12_patches_bi,r12_patches_bi_internal,0);
                        
                        %%% match the second connection
                        else

                            %%% look for angle over second connection
                            theta2_init = 0;
                            for ai = 1:o.nangle
                                if o.angles_bis(2,ai) == b0 && o.angles_ibs(2,ai) == ib_conn2_b0
                                    if o.angles_bis(3,ai) == bi && o.angles_ibs(3,ai) == ib_conn2_b1
                                        r12_patches_b0_angle2 = o.bs(b0).r(:,o.angles_ibs(2,ai)) - o.bs(b0).r(:,o.angles_ibs(1,ai));
                                        theta2_init = o.angles_theta_init(ai);
                                        break
                                    end
                                end
                                if o.angles_bis(3,ai) == b0 && o.angles_ibs(3,ai) == ib_conn_b0
                                    if o.angles_bis(2,ai) == bi && o.angles_ibs(2,ai) == ib_conn_b1
                                        r12_patches_b0_angle2 = o.bs(b0).r(:,o.angles_ibs(4,ai)) - o.bs(b0).r(:,o.angles_ibs(3,ai));
                                        theta2_init = o.angles_theta_init(ai);
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
                                r12_conns_b0 = ars.unitVector(o.bs(b0).r(:,ib_conn2_b0) - o.bs(b0).r(:,ib_conn_b0));
                                r12_perp = (2*randi(2)-3)*cross(ars.unitVector(r12_patches_b0),r12_conns_b0);
                                r12_conn = r12_conn_mag*( cosd(180-theta_init)*ars.unitVector(r12_patches_b0) + sind(180-theta_init)*r12_perp );

                                %%% block direction
                                r12_perp = (2*randi(2)-3)*cross(ars.unitVector(r12_conn),r12_conns_b0);
                                r12_patches_bi = cosd(180-theta_init)*ars.unitVector(r12_conn) + sind(180-theta_init)*r12_perp;
                                r12_patches_bi_internal = ars.unitVector(o.bs(bi).r_internal(:,ib_end_b1) - o.bs(bi).r_internal(:,ib_conn_b1));
                                [R,failed_rot,r12_conn2] = o.optimizeR(r12_patches_bi,r12_patches_bi_internal,b0,bi,ib_conn_b0,ib_conn_b1,ib_conn2_b0,ib_conn2_b1,r12_conn,r12_eq_conn2,tolerance_rot);
                                if failed_rot
                                    attempts_rot = attempts_rot + 1;
                                    continue
                                end

                                %%% check second angle
                                if theta2_init ~= 0
                                    if abs(dot(ars.unitVector(r12_patches_b0_angle2), ars.unitVector(r12_conn2)) - cosd(180-theta2_init)) > tolerance_rot
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
                    r_source = o.bs(b0).r(:,ib_conn_b0);
                    [o.bs(bi),failed_block,r_other_origami] = o.bs(bi).init_positions(p,r_source,r12_conn,ib_conn_b1,R,r_other_origami);
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
                if are_conns_overstretched(o,p,U_overstretched)
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
        function [R,failed,r12_conn2] = optimizeR(o,r12_patches_bi,r12_patches_bi_internal,b0,bi,ib_conn_b0,ib_conn_b1,ib_conn2_b0,ib_conn2_b1,r12_conn,r12_eq_conn2,tolerance)
            max_iterations = 1000;
            iteration = 0;
            initializing = true;
            finding_min = true;
            phi = rand*360;
            phi_prev = phi;
            step = 10;
            while true
                if iteration == max_iterations
                    failed = true;
                    break
                end
                iteration = iteration + 1;
                R = o.getRotation(r12_patches_bi,r12_patches_bi_internal,phi);
                r_conn_bi = o.bs(b0).r(:,ib_conn_b0) + r12_conn;
                r_conn_bi_internal = o.bs(bi).r_internal(:,ib_conn_b1);
                r_conn2_bi_internal = o.bs(bi).r_internal(:,ib_conn2_b1);
                r_conn2_bi = r_conn_bi + R'*(r_conn2_bi_internal-r_conn_bi_internal);
                r_conn2_b0 = o.bs(b0).r(:,ib_conn2_b0);
                r12_conn2 = r_conn2_bi-r_conn2_b0;
                d = norm(r12_conn2);
                if initializing == true
                    if phi - phi_prev == step
                        m = (d-d_prev)/step;
                        if m > 0
                            step = -step;
                        end
                        m_prev = m;
                        initializing = false;
                    end
                elseif abs(d-r12_eq_conn2) < tolerance
                    failed = false;
                    break
                elseif finding_min == true
                    if d-r12_eq_conn2 < 0
                        finding_min = false;
                        d = d_prev;
                        phi = phi - step;
                    else
                        m = (d-d_prev)/step;
                        if sign(m) ~= sign(m_prev)
                            step = -step/2;
                        end
                        m_prev = m;
                    end
                else
                    if sign(d-r12_eq_conn2) ~= sign(d_prev-r12_eq_conn2)
                        step = -step/2;
                    end
                end
                d_prev = d;
                phi = phi + step;
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
            R0 = B*A';
            u = b; U = [0 -u(3) u(2); u(3) 0 -u(1); -u(2) u(1) 0];
            Rspin = eye(3) + sind(phi)*U + (1-cosd(phi))*(U*U);
            R = Rspin * R0;
        end
    end
end