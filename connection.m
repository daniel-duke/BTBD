%%% connection class for BTBD 
classdef connection
    properties
        bis          % block indices
        is           % indices of connected beads in block
        r12_eq       % equilibrium separation
    end

    methods
        %%% constructor
        function c = connection(bi1,i1,bi2,i2,r12_eq)
            if nargin > 0
                c.bis = [bi1,bi2];
                c.is = [i1,i2];
                c.r12_eq = r12_eq;
            end
        end

    end
end