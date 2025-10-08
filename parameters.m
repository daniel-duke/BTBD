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
        U_overstretched
    end

    methods
        %%% constructor
        function p = parameters(nstep_eq,shrink_ratio,nstep_prod,dump_every,...
                                dt,dbox,verlet_skin,neigh_every,react_every,...
                                T,sigma,epsilon,r12_helix,r12_bead,U_overstretched)
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
                p.U_overstretched = 6.96*U_overstretched;
            end
        end
    end
end