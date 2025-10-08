%%% potential class for BTBD
classdef angle_pot
    properties
        theta_eq        % equilibrium angle
        k_theta         % spring constant
        index           % angle type
    end

    methods
        %%% constructor
        function apot = angle_pot(theta_eq,k_theta,index)
            if nargin > 0
                apot.theta_eq = theta_eq;
                apot.k_theta = 6.96*k_theta;
                apot.index = index;
            end
        end

    end
end