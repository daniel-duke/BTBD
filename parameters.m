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
        r12_helix
        r12_bead
        k_x_conn
        k_x_linker
        k_theta
    end

    methods
        %%% constructor
        function p = parameters(nstep_eq,shrink_ratio,nstep_prod,dump_every,...
                                dt,dbox,verlet_skin,neigh_every,react_every,...
                                T,sigma,epsilon,r12_helix,r12_bead,...
                                k_x_conn,k_x_linker,k_theta)
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
                p.r12_helix = r12_helix;
                p.r12_bead = r12_bead;
                p.k_x_conn = 6.96*k_x_conn;
                p.k_x_linker = 6.96*k_x_linker;
                p.k_theta = k_theta;
            end
        end
    end
end