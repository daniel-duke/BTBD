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


        %%% check connections
        function overstretched = check_connections(o)
            overstretched = false;
            for ci = 1:length(o.cs)
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
                if norm(r1-r2) > 2*o.cs(ci).r12_eq
                    overstretched = true;
                    return
                end
            end
        end
        

        %%% initialize origami positions
        function [o,overlap] = init(o,p,r_other)
            max_attempts = 100;
            box_ratio = 3/4;

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

            if o.design == "1B"
                [o,overlap] = init_1B(o,p,r_other,max_attempts,box_ratio);
            
            elseif o.design == "2B"
                [o,overlap] = init_2B(o,p,r_other,max_attempts,box_ratio);

            elseif o.design == "3B"
                [o,overlap] = init_3B(o,p,r_other,max_attempts,box_ratio);

            else
                error("Unknown origami design.")

            end
        end


        %%% update overall origami positions
        function r = update_r(o)
            r = [];
            if ~isempty(o.ts); r = [r o.ts.r]; end
            if ~isempty(o.bs); r = [r o.bs.r]; end
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Initialization Functions %%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%% initialize positions for 1 block design
        function [o,overlap] = init_1B(o,p,r_other,max_attempts,box_ratio)
    
            %%% initialize positions loop (success only if no overlap)
            r_other_initial = r_other;
            for attempts = 1:max_attempts

                %%% place block
                r_conn = (rand(3,1)-0.5)*p.dbox*box_ratio;
                [o.bs,overlap] = o.bs.init(p,1,r_conn,box_muller,r_other);
                if overlap == true
                    r_other = r_other_initial;
                    continue
                end

                %%% making it this far means success
                break
            end
        end


        %%% initialize positions for 2 block
        function [o,overlap] = init_2B(o,p,r_other,max_attempts,box_ratio)

            %%% initialize positions loop (success only if no overlap)
            r_other_initial = r_other;
            for attempts = 1:max_attempts

                %%% place first block
                bi = 1;
                r12_b1 = box_muller;
                r_conn = (rand(3,1)-0.5)*p.dbox*box_ratio;
                i_conn = 1;
                [o.bs(bi),overlap] = o.bs(1).init(p,i_conn,r_conn,r12_b1,r_other);
                if overlap == true
                    r_other = r_other_initial;
                    continue
                end
                r_other = ars.my_horzcat(r_other,o.bs(bi).r(:,o.bs(bi).n_r));

                %%% place second block
                bi = 2;
                r12_b2 = cross(r12_b1,box_muller);
                r_conn = r_conn + p.sigma/sqrt(2)*ars.unit_vector(-r12_b1) + p.sigma/sqrt(2)*ars.unit_vector(r12_b2);
                i_conn = 1;
                [o.bs(bi),overlap] = o.bs(bi).init(p,i_conn,r_conn,r12_b2,r_other);
                if overlap == true
                    r_other = r_other_initial;
                    continue
                end

                %%% check connections
                overstretched = o.check_connections;
                if overstretched == true
                    overlap = true;
                    r_other = r_other_initial;
                    continue
                end

                %%% making it this far means success
                break
            end
        end


        %%% initialize positions for 3 block design
        function [o,overlap] = init_3B(o,p,r_other,max_attempts,box_ratio)

            %%% initialize positions loop (success only if no overlap)
            r_other_initial = r_other;
            for attempts = 1:max_attempts

                %%% place first block
                bi = 1;
                r12 = box_muller;
                r_conn = (rand(3,1)-0.5)*p.dbox*box_ratio;
                i_conn = o.bs(bi).interpret_loc(o.conn{1});
                [o.bs(bi),overlap] = o.bs(bi).init(p,i_conn,r_conn,r12,r_other);
                if overlap == true
                    r_other = r_other_initial;
                    continue
                end
                r_other = ars.my_horzcat(r_other,o.bs(bi).r(:,o.bs(bi).n_r));

                %%% place second block
                bi = 2;
                r12 = -r12;
                r_conn = r_conn + p.r12_eq_ghost.*ars.unit_vector(r12);
                i_conn = o.bs(bi).interpret_loc(o.conn{1});
                [o.bs(bi),overlap] = o.bs(bi).init(p,i_conn,r_conn,r12,r_other);
                if overlap == true
                    r_other = r_other_initial;
                    continue
                end
                r_other = ars.my_horzcat(r_other,o.bs(bi).r(:,o.bs(bi).n_r));

                %%% place third block
                bi = 3;
                r12 = cross(r12,box_muller);
                r_conn = r_conn + p.r12_eq_ghost.*ars.unit_vector(r12);
                i_conn = o.bs(bi).interpret_loc(o.conn{1});
                [o.bs(bi),overlap] = o.bs(bi).init(p,i_conn,r_conn,r12,r_other);
                if overlap == true
                    r_other = r_other_initial;
                    continue
                end

                %%% check connections
                overstretched = check_connections(o,p);
                if overstretched == true
                    overlap = true;
                    r_other = r_other_initial;
                    continue
                end

                %%% making it this far means success
                break
            end
        end
        
    end
end