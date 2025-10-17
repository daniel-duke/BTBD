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
        comm_cutoff
    end

    methods
        %%% constructor
        function p = parameters(p_input)
            if nargin > 0
                p.nstep = p_input.nstep;
                p.nstep_relax = p_input.nstep_relax;
                p.dump_every = p_input.dump_every;
                p.dbox = p_input.dbox;
                p.shrink_ratio = p_input.shrink_ratio;
                p.dt = p_input.dt;
                p.verlet_skin = p_input.verlet_skin;
                p.neigh_every = p_input.neigh_every;
                p.react_every = p_input.react_every;
                p.T = p_input.T;
                p.sigma = p_input.sigma;
                p.epsilon = 6.96*p_input.epsilon;
                p.r12_cut_WCA = p.sigma*2^(1/6);
                p.r12_helix = p_input.r12_helix;
                p.r12_bead = p_input.r12_bead;
                p.U_strained = 6.96*p_input.U_strained;
                p.comm_cutoff = p_input.comm_cutoff;
            end
        end
    end
end