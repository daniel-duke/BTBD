%%% parameters class for BTBD
classdef parameters
    properties
        nstep
        nstep_relax
        dump_every
        dbox
        shrink_ratio
        dt
        verlet_skin
        neigh_every
        react_every
        T
        sigma
        epsilon
        r12_cut_WCA
        r12_helix
        r12_bead
        U_strained
    end

    methods
        %%% constructor
        function p = parameters(p_in)
            if nargin > 0
                p.nstep = p_in.nstep;
                p.nstep_relax = p_in.nstep_relax;
                p.dump_every = p_in.dump_every;
                p.dbox = p_in.dbox;
                p.shrink_ratio = p_in.shrink_ratio;
                p.dt = p_in.dt;
                p.verlet_skin = p_in.verlet_skin;
                p.neigh_every = p_in.neigh_every;
                p.react_every = p_in.react_every;
                p.T = p_in.T;
                p.sigma = p_in.sigma;
                p.epsilon = 6.96*p_in.epsilon;
                p.r12_cut_WCA = p.sigma*2^(1/6);
                p.r12_helix = p_in.r12_helix;
                p.r12_bead = p_in.r12_bead;
                p.U_strained = 6.96*p_in.U_strained;
            end
        end
    end
end