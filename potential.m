%%% potential class for BTBD
classdef potential
    properties
        style           % bond style
        r12_eq          % equilibrium separation
        params          % other bond parameters
        index           % bond type
    end

    methods
        %%% constructor
        function pot = potential(style,r12_eq,params,index)
            if nargin > 0
                pot.style = style;
                pot.r12_eq = r12_eq;
                pot.params = params;
                pot.index = index;
                pot = pot.adjust_units();
            end
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
        
        %%% convert energy parameters to kcal/mol
        function pot = adjust_units(pot)
            
            %%% harmonic bond
            if pot.style == "harmonic"
                pot.params(1) = 6.96*pot.params(1);

            %%% cut harmonic bond
            elseif pot.style == "harmonic/shift/cut"
                pot.params(1) = 6.96*pot.params(1);
                pot.params(2) = 6.96*pot.params(2);

            end
        end


        %%% calculate separation distance at energy
        function r12_mag = calc_separation(pot,U)

            %%% dummy bond
            if pot.style == "zero"
                r12_mag = pot.r12_eq;
        
            %%% harmonic bond
            elseif pot.style == "harmonic"
                k_x = pot.params(1);
                r12_mag = pot.r12_eq + sqrt(2*U/k_x);
        
            %%% cut harmonic bond
            elseif pot.style == "harmonic/shift/cut"
                k_x = pot.params(1);
                U_cut = pot.params(2);
                if U < U_cut
                    r12_mag = pot.r12_eq + sqrt(2*U/k_x);
                else
                    r12_mag = pot.r12_eq + sqrt(2*U_cut/k_x);
                end

            %%% error
            else
                error("Unrecognized potential.")
            end
        end

        %%% calculate energy at separation distance
        function U = calc_energy(pot,r12_mag)

            %%% dummy bond
            if pot.style == "zero"
                U = 0;
        
            %%% harmonic bond
            elseif pot.style == "harmonic"
                k_x = pot.params(1);
                U = k_x/2*(r12_mag-pot.r12_eq)^2;
        
            %%% cut harmonic bond
            elseif pot.style == "harmonic/shift/cut"
                k_x = pot.params(1);
                U_cut = pot.params(2);
                r12_cut = pot.r12_eq + sqrt(2*U_cut/k_x);
                if r12_mag < r12_cut
                    U = k_x/2*(r12_mag-pot.r12_eq)^2;
                else
                    U = U_cut;
                end
        
            %%% error
            else
                error("Unrecognized potential.")
            end
        end
        
        
        %%% write coeffecients for potential
        function write_potential(pot,f)
            fprintf(f,strcat(...
                "bond_coeff      ",...
                num2str(pot.index), " ",...
                pot.style));

            %%% dummy bond
            if pot.style == "zero"
                fprintf(f,"\n");
        
            %%% harmonic bond
            elseif pot.style == "harmonic"
                k_x = pot.params(1);
                fprintf(f,strcat(" ",...
                    num2str(k_x/2)," ",...
                    num2str(pot.r12_eq),"\n"));
        
            %%% cut harmonic bond
            elseif pot.style == "harmonic/shift/cut"
                k_x = pot.params(1);
                U_cut = pot.params(2);
                r12_cut = pot.r12_eq + sqrt(2*U_cut/k_x);
                fprintf(f,strcat(" ",...
                    num2str(U_cut)," ",...
                    num2str(pot.r12_eq)," ",...
                    num2str(r12_cut)," ","\n"));
        
            %%% error
            else
                error("Unrecognized potential.")
            end
        end

    end
end