%%% parameters class for BTBD
classdef parameters
    properties
        nstep_eq
        shrink_ratio
        nstep_prod
        dump_every
        dt
        dbox
        verlet_skin
        neigh_every
        react_every
        T
        sigma
        epsilon
        r12_cut_WCA
        r12_eq_block
        r12_adj_block
        k_x_conn
        r12_eq_linker
        k_x_linker
    end

    methods
        %%% constructor
        function p = parameters(nstep_eq,shrink_ratio,nstep_prod,dump_every,...
                                dt,dbox,verlet_skin,neigh_every,react_every,...
                                T,sigma,epsilon,r12_eq_block,r12_adj_block,...
                                k_x_conn,r12_eq_linker,k_x_linker)
            if nargin > 0
                p.nstep_eq = nstep_eq;
                p.shrink_ratio = shrink_ratio;
                p.nstep_prod = nstep_prod;
                p.dump_every = dump_every;
                p.dt = dt;
                p.dbox = dbox;
                p.verlet_skin = verlet_skin;
                p.neigh_every = neigh_every;
                p.react_every = react_every;
                p.T = T;
                p.sigma = sigma;
                p.epsilon = 6.96*epsilon;
                p.r12_cut_WCA = p.sigma*2^(1/6);
                p.r12_eq_block = r12_eq_block;
                p.r12_adj_block = r12_adj_block;
                p.k_x_conn = 6.96*k_x_conn;
                p.r12_eq_linker = r12_eq_linker;
                p.k_x_linker = 6.96*k_x_linker;
            end
        end
    end
end