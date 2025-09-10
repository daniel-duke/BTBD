%%% origami class for BTBD
classdef origami
    properties
        bs          % block objects
        cs          % internal connections
        conn        % internal connection info
        n           % total number of particles
        r           % total positions
        nconn       % number of connections
        bi_conn     % connection block indices
        i_conn      % connection bead indices
        get_bi      % map pi to bi (block index)
        get_i       % map pi to i (bead index)
        get_pi      % map bi and i to pi (index within origami)
    end

    methods
        %%% constructor
        function o = origami()
            if nargin > 0
                o.bs = [];
                o.cs = [];
                o.n = sum([bs.n]);
                o.r = zeros(3,o.n);
                o.nconn = 0;
                o.bi_conn = zeros(2,0);
                o.i_conn = zeros(2,0);
                [o.get_bi,o.get_i,o.get_pi] = map_indices(o);
            end
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Utility Functions %%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%% add connection to origami
        function o = add_conn(o,bi1,loc1,bi2,loc2)
            o.nconn = o.nconn + 1;
            o.bi_conn = [ o.bi_conn [bi1;bi2] ];
            i1 = o.bs(bi1).interpret_loc(loc1);
            i2 = o.bs(bi2).interpret_loc(loc2);
            o.i_conn = [ o.i_conn [i1;i2] ];
        end


        %%% map between particle indices and block/bead indices
        function [get_bi,get_i,get_pi] = map_indices(o)
            get_bi = zeros(1,o.n);
            get_i = zeros(1,o.n);
            get_pi = zeros(length([o.bs]),o.bs.n);
            pi = 0;
            for bi = 1:length([o.bs])
                for i = 1:o.bs(bi).n
                    pi = pi+1;
                    get_bi(pi) = bi;
                    get_i(pi) = i;
                    get_pi(bi,i) = pi;
                end
            end
        end

        
        %%% check if particle is real
        function result = is_real(o,pi)
            if o.get_i(pi) <= o.bs(o.get_bi(pi)).n_r
                result = true;
            else
                result = false;
            end
        end


        %%% update overall origami positions
        function r = update_r(o)
            r = [o.bs.r];
        end


        %%% ensure connections are not too strained
        function overstretch = check_connections(o,p,U_overstretched)
            for ci = 1:length(o.cs)
                r12_eq = o.cs(ci).r12_eq;
                r12_overstretched = r12_eq + sqrt(2*U_overstretched/p.k_x_conn);
                r1 = o.bs(o.cs(ci).bis(1)).r(:,o.cs(ci).is(1));
                r2 = o.bs(o.cs(ci).bis(2)).r(:,o.cs(ci).is(2));
                if norm(r1-r2) > r12_overstretched
                    overstretch = true;
                    return
                end
            end
            overstretch = false;
        end
        

        %%% initialize origami positions
        function [o,fail,r_other] = init(o,p,r_other)
            max_attempts = 1000;
            box_ratio = 3/4;
            U_overstretched = 100;

            %%% create connections
            o.cs = connection.empty;
            nconn = length(o.conn);
            for i = 1:nconn
                bi1 = o.conn{i}{1};
                bi2 = o.conn{i}{3};
                o.cs(i) = connection(...
                    bi1,o.bs(bi1).interpret_loc(o.conn{i}{2}),...
                    bi2,o.bs(bi2).interpret_loc(o.conn{i}{4}),...
                    o.conn{i}{5});
            end

            %%% store positions
            r_other_initial = r_other;

            %%% origami attempt loop
            attempts = 0;
            while true

                %%% check origami attempts
                if attempts == max_attempts
                    fail = true;
                    break
                end

                %%% add block
                r_source = (rand(3,1)-0.5)*p.dbox*box_ratio;
                [o.bs(1),overlap,r_other] = o.bs(1).init(p,1,r_source,0,r_other);

                %%% reset if overlap
                if overlap == true
                    attempts = attempts + 1;
                    r_other = r_other_initial;
                    continue
                end

                %%% loop over remaining blocks
                for bi = 2:length(o.bs)

                    %%% find connection between block and previous block
                    found_connection = false;
                    for ci = 1:length(o.cs)
                        if o.cs(ci).bis(1) == bi-1 && o.cs(ci).bis(2) == bi
                            i_conn_prev = o.cs(ci).is(1);
                            i_conn_curr = o.cs(ci).is(2);
                            r12_eq_conn = o.cs(ci).r12_eq;
                            found_connection = true;
                            break
                        end
                        if o.cs(ci).bis(2) == bi-1 && o.cs(ci).bis(1) == bi
                            i_conn_prev = o.cs(ci).is(2);
                            i_conn_curr = o.cs(ci).is(1);
                            r12_eq_conn = o.cs(ci).r12_eq;
                            found_connection = true;
                            break
                        end
                    end
                    if found_connection == false
                        error("Cound not find connection between block and previous block.")
                    end

                    %%% add block
                    r_source = o.bs(bi-1).r(:,i_conn_prev);
                    [o.bs(bi),overlap,r_other] = o.bs(bi).init(p,i_conn_curr,r_source,r12_eq_conn,r_other);
                    if overlap == true
                        break
                    end
                end

                %%% check connections
                overstretch = check_connections(o,p,U_overstretched);

                %%% reset if block overlap or connection overstretch
                if overlap == true || overstretch == true
                    attempts = attempts + 1;
                    r_other = r_other_initial;
                    continue
                end

                %%% origami succesfully initialized
                fail = false;
                break
            end
        end
        
    end
end