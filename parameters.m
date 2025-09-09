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
        kBT
        sigma
        epsilon
        r12_cut_WCA
        r12_eq_block
        r12_adj_block
        r12_eq_tether
        k_x_tether
        k_theta
        k_x_conn
        r12_eq_linker
        k_x_linker
        r12_cut_linker
    end

    methods
        %%% constructor
        function p = parameters(nstep_eq,shrink_ratio,nstep_prod,dump_every,...
                                dt,dbox,verlet_skin,neigh_every,react_every,...
                                T,sigma,epsilon,r12_eq_block,r12_adj_block,r12_eq_tether,k_x_tether,Lp_tether,...
                                k_x_conn,r12_eq_linker,k_x_linker,U_cut_linker)
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
                p.kBT = 0.0138*p.T;
                p.sigma = sigma;
                p.epsilon = 6.96*epsilon;
                p.r12_cut_WCA = p.sigma*2^(1/6);
                p.r12_eq_block = r12_eq_block;
                p.r12_adj_block = r12_adj_block;
                p.r12_eq_tether = r12_eq_tether;
                p.k_x_tether = 6.96*k_x_tether;
                p.k_theta = Lp_tether*0.0138/p.r12_eq_tether;
                p.k_x_conn = 6.96*k_x_conn;
                p.r12_eq_linker = r12_eq_linker;
                p.k_x_linker = 6.96*k_x_linker;
                p.r12_cut_linker = p.r12_eq_linker + sqrt(2*U_cut_linker*6.96/p.k_x_linker);
            end
        end
    end
end