%%% origami class for BTBD
classdef origami
    properties
        design      %structure type
        ts          %tether objects
        bs          %block objects
        cs          %internal connections
        conn        %internal connection info
        n           %total number of particles
        r           %total positions
        get_dom     %map pi to domain identity (1 for tether, 2 for block, 3 for patch)
        get_domi    %map pi to domain index (which tether, or which block)
        get_i       %map pi to bead index (index within domain)
        get_pi      %map dom, domi, and i to pi (index within origami)
    end

    methods
        %%% constructor
        function o = origami(design,ts,bs,conn)
            if nargin > 0
                o.design = design;
                o.ts = ts;
                o.bs = bs;
                o.cs = [];
                o.conn = conn;
                o.n = sum([ts.n]) + sum([bs.n]);
                o.r = zeros(3,o.n);
                [o.get_dom,o.get_domi,o.get_i,o.get_pi] = map_indices(o);
            end
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Utility Functions %%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%% update overall origami positions
        function r = update_r(o)
            r = [];
            if ~isempty(o.ts); r = [r o.ts.r]; end
            if ~isempty(o.bs); r = [r o.bs.r]; end
        end


        %%% map between particle indices and domain/bead indices
        function [get_dom,get_domi,get_i,get_pi] = map_indices(o)
            max_domi = max([length([o.ts]),length([o.bs])]);
            if ~isempty(o.ts); max_i = max([o.ts.n]); else; max_i = 0; end
            if ~isempty(o.bs); max_i = max([max_i,[o.bs.n]]); end
            get_dom = zeros(1,o.n);
            get_domi = zeros(1,o.n);
            get_i = zeros(1,o.n);
            get_pi = zeros(2,max_domi,max_i);
            pi = 0;
            for ti = 1:length([o.ts])
                for i = 1:o.ts(ti).n
                    pi = pi+1;
                    get_dom(pi) = 1;
                    get_domi(pi) = ti;
                    get_i(pi) = i;
                    get_pi(1,ti,i) = pi;
                end
            end
            for bi = 1:length([o.bs])
                for i = 1:o.bs(bi).n
                    pi = pi+1;
                    get_dom(pi) = 2;
                    get_domi(pi) = bi;
                    get_i(pi) = i;
                    get_pi(2,bi,i) = pi;
                end
            end
        end


        %%% ensure connections are not too strained
        function overstretch = check_connections(o,p,U_overstretched)
            for ci = 1:length(o.cs)
                r12_eq = o.cs(ci).r12_eq;
                r12_overstretched = r12_eq + sqrt(2*U_overstretched/p.k_x_ghost);
                if o.cs(ci).doms(1) == 1
                    r1 = o.ts(o.cs(ci).domis(1)).r(:,o.cs(ci).is(1));
                else
                    r1 = o.bs(o.cs(ci).domis(1)).r(:,o.cs(ci).is(1));
                end
                if o.cs(ci).doms(2) == 1
                    r2 = o.ts(o.cs(ci).domis(2)).r(:,o.cs(ci).is(2));
                else
                    r2 = o.bs(o.cs(ci).domis(2)).r(:,o.cs(ci).is(2));
                end
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
            U_overstretched = 10;

            %%% create connections
            o.cs = connection.empty;
            nconn = length(o.conn);
            for i = 1:nconn
                bi1 = o.conn{i}{1};
                bi2 = o.conn{i}{3};
                o.cs(i) = connection(...
                    2,bi1,o.bs(bi1).interpret_loc(o.conn{i}{2}),...
                    2,bi2,o.bs(bi2).interpret_loc(o.conn{i}{4}),...
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
                        if o.cs(ci).domis(1) == bi-1 && o.cs(ci).domis(2) == bi
                            i_conn_prev = o.cs(ci).is(1);
                            i_conn_curr = o.cs(ci).is(2);
                            r12_eq_conn = o.cs(ci).r12_eq;
                            found_connection = true;
                            break
                        end
                        if o.cs(ci).domis(2) == bi-1 && o.cs(ci).domis(1) == bi
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