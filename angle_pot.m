%%% angle potential class for BTBD
classdef angle_pot
    properties
        label           % name
        style           % angle style
        theta_eq        % equilibrium angle
        params          % parameters object
        index           % angle type (numeric)
    end

    methods
        %%% constructor
        function apot = angle_pot(style,theta_eq,params_input,index)
            if nargin > 0
                apot.label = style;
                apot.theta_eq = theta_eq;
                apot.index = index;
                apot = apot.set_params(params_input);
            end
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%% put names to parameters, convert energy to kcal/mol
        function apot = set_params(apot,params_input)

            %%% harmonic angle
            if apot.label == "harmonic"
                apot.style = "harmonic";
                apot.params.k_theta = 6.96*params_input(1);

            %%% cosine angle
            elseif apot.label == "cosine"
                apot.style = "cosine/delta";
                apot.params.k_theta = 6.96*params_input(1);

            %%% error
            else
                error("Unknown angle type.")
            end
        end


        %%% calculate energy at separation distance
        function U = calc_energy(apot,theta,theta_ref)
            arguments
                apot; theta
                theta_ref = apot.theta_eq
            end

            %%% harmonic angle
            if apot.label == "harmonic"
                U = apot.params.k_theta/2*(theta-theta_ref)^2;

            %%% cosine angle
            elseif apot.label == "cosine"
                U = apot.params.k_theta*(1-cosd(theta-theta_ref));

            end
        end
        
        
        %%% write coeffecients for potential
        function write_potential(apot,f,is_hybrid)
            fprintf(f,strcat(...
                "angle_coeff     ",...
                num2str(apot.index)));
            if is_hybrid
                fprintf(f," " + apot.style);
            end

            %%% harmonic angle
            if apot.label == "harmonic"
                fprintf(f,strcat(" ",...
                    num2str(apot.params.k_theta/2)," ",...
                    num2str(apot.theta_eq),"\n"));

            %%% harmonic angle
            elseif apot.label == "cosine"
                fprintf(f,strcat(" ",...
                    num2str(apot.params.k_theta)," ",...
                    num2str(apot.theta_eq),"\n"));
                
            end
        end

    end
end