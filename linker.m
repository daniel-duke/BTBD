%%% linker class for BTBD
classdef linker
    properties
        pot_index       % bond type
        r12_cut         % cutoff radius
        apot_index      % angle type
        theta_min       % minimum theta
        theta_max       % maximum theta
        index           % linker type
    end

    methods
        %%% constructor
        function l = linker(pot_index,r12_cut,index)
            if nargin > 0
                l.pot_index = pot_index;
                l.r12_cut = r12_cut;
                l.apot_index = 0;
                l.theta_min = 0;
                l.theta_max = 180;
                l.index = index;
            end
        end

    end
end