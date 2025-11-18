%%% parameters class for BTBD
classdef parameters
    properties
        rseed
        rseed_lmp
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
        r12_xy
        r12_z
        mass
        U_strained
        nABADtype
        n_vis
    end

    methods
        %%% constructor
        function p = parameters(p_input)
            if nargin > 0
                p.rseed = p_input.rseed;
                p.rseed_lmp = p_input.rseed_lmp;
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
                p.r12_xy = p_input.r12_xy;
                p.r12_z = p_input.r12_z;
                p.mass = p_input.mass;
                p.U_strained = 6.96*p_input.U_strained;
                p.nABADtype = zeros(1,4);
                p.n_vis = 0;
            end
        end

    end
    properties (Constant)

        %%% default values for input file parameters
        defaults = struct(...
            'rseed',        NaN, ...        % none          - random seed for BTBD
            'rseed_lmp',    1,   ...        % none          - random seed for LAMMPS
            'nstep',        NaN, ...        % steps         - number of production steps
            'nstep_relax',  1e5, ...        % steps         - number of relaxation/shrink steps
            'dump_every',   1e4, ...        % steps         - number of steps between dumps
            'dbox',         NaN, ...        % nm            - periodic boundary diameter
            'shrink_ratio', 1,   ...        % none          - box compression (final/initial)
            'dt',           NaN, ...        % ns            - integration time step
            'verlet_skin',  4,   ...        % nm            - width of neighbor list skin
            'neigh_every',  1e1, ...        % steps         - how often to consider updating neighbor list
            'react_every',  1e1, ...        % steps         - how often to check for linker hybridization
            'T',            300, ...        % K             - temperature
            'sigma',        NaN, ...        % nm            - WCA distance parameter
            'epsilon',      NaN, ...        % kcal/mol      - WCA energy parameter
            'r12_xy',       NaN, ...        % nm            - structural bead separation in block xy-plane
            'r12_z',        NaN, ...        % nm            - structural bead separation along block z-axis
            'mass',         1,   ...        % none          - mass of structural beads
            'U_strained',   10   ...        % kcal/mol      - max energy for initialized bonds and angles
        );

        %%% old and new parameter names
        name_changes = struct(...
            'r12_helix', 'r12_xy',...
            'r12_bead', 'r12_z'...
        );

    end
end