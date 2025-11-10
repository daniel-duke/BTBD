%%% angle potential class for BTBD
classdef angle_pot
    properties
        style           % angle style (BTBD notation)
        style_lmp       % angle style (LAMMPS notation)
        theta_eq        % equilibrium angle
        params          % parameters object
        index           % angle type (numeric)
    end

    methods
        %%% constructor
        function apot = angle_pot(style,theta_eq,params_input,index)
            if nargin > 0
                apot.style = style;
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
            if apot.style == "harmonic"
                apot.style_lmp = "harmonic";
                apot.params.k_theta = 6.96*params_input(1);

            %%% cosine angle
            elseif apot.style == "cosine"
                apot.style_lmp = "cosine/delta";
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
            if apot.style == "harmonic"
                U = apot.params.k_theta/2*(theta-theta_ref)^2;

            %%% cosine angle
            elseif apot.style == "cosine"
                U = apot.params.k_theta*(1-cosd(theta-theta_ref));

            end
        end
        
        
        %%% write coeffecients for potential
        function write_potential(apot,f,is_hybrid)
            fprintf(f,strcat(...
                "angle_coeff     ",...
                num2str(apot.index)));
            if is_hybrid
                fprintf(f," " + apot.style_lmp);
            end

            %%% harmonic angle
            if apot.style == "harmonic"
                fprintf(f,strcat(" ",...
                    ars.fstring(apot.params.k_theta/2,0,2)," ",...
                    ars.fstring(apot.theta_eq,0,2),"\n"));

            %%% harmonic angle
            elseif apot.style == "cosine"
                fprintf(f,strcat(" ",...
                    ars.fstring(apot.params.k_theta,0,2)," ",...
                    ars.fstring(apot.theta_eq,0,2),"\n"));
                
            end
        end

    end
end