%%% connection class for BTBD 
classdef connection
    properties
        doms         % domain types (1 for tether, 2 for block)
        domis        % indices of domain
        is           % indices of connected beads in domin
        r12_eq       % equilibrium separation
    end

    methods
        %%% constructor
        function c = connection(dom1,domi1,i1,dom2,domi2,i2,r12_eq)
            if nargin > 0
                c.doms = [dom1,dom2];
                c.domis = [domi1,domi2];
                c.is = [i1,i2];
                c.r12_eq = r12_eq;
            end
        end

    end
end