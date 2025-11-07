%%% bond potential class for BTBD
classdef bond_pot
    properties
        style           % bond style (BTBD notation)
        style_lmp       % bond style (LAMMPS notation)
        r12_eq          % equilibrium separation
        params          % parameters object
        index           % bond type (numeric)
    end

    methods
        %%% constructor
        function pot = bond_pot(style,r12_eq,params_input,index)
            if nargin > 0
                pot.style = style;
                pot.r12_eq = r12_eq;
                pot.index = index;
                pot = pot.set_params(params_input);
            end
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%% put names to parameters, convert energy to kcal/mol
        function pot = set_params(pot,params_input)
            
            %%% dummy bond
            if pot.style == "zero"
                pot.style_lmp = "zero";

            %%% harmonic bond
            elseif pot.style == "harmonic"
                pot.style_lmp = "harmonic";
                pot.params.k_x = 6.96*params_input(1);

            %%% cut harmonic bond
            elseif pot.style == "harmonic_cut"
                pot.style_lmp = "harmonic/shift/cut";
                pot.params.k_x = 6.96*params_input(1);
                pot.params.U_cut = 6.96*params_input(2);
                pot.params.r12_cut = pot.r12_eq + sqrt(2*pot.params.U_cut/pot.params.k_x);

            %%% error
            else
                error("Unknown bond type.")
            end

        end


        %%% calculate separation distance at energy
        function r12_mag = calc_separation(pot,U)

            %%% dummy bond
            if pot.style == "zero"
                r12_mag = pot.r12_eq;
        
            %%% harmonic bond
            elseif pot.style == "harmonic"
                r12_mag = pot.r12_eq + sqrt(2*U/pot.params.k_x);
        
            %%% cut harmonic bond
            elseif pot.style == "harmonic_cut"
                if U < pot.params.U_cut
                    r12_mag = pot.r12_eq + sqrt(2*U/pot.params.k_x);
                else
                    r12_mag = pot.r12_eq + sqrt(2*pot.params.U_cut/pot.params.k_x);
                end

            end
        end


        %%% calculate energy at separation distance
        function U = calc_energy(pot,r12_mag)

            %%% dummy bond
            if pot.style == "zero"
                U = 0;
        
            %%% harmonic bond
            elseif pot.style == "harmonic"
                U = pot.params.k_x/2*(r12_mag-pot.r12_eq)^2;
        
            %%% cut harmonic bond
            elseif pot.style == "harmonic_cut"
                if r12_mag < r12_cut
                    U = pot.params.k_x/2*(r12_mag-pot.r12_eq)^2;
                else
                    U = pot.params.U_cut;
                end

            end
        end
        
        
        %%% write coeffecients for potential
        function write_potential(pot,f,is_hybrid)
            fprintf(f,strcat(...
                "bond_coeff      ",...
                num2str(pot.index)));
            if is_hybrid
                fprintf(f," " + pot.style_lmp);
            end

            %%% dummy bond
            if pot.style == "zero"
                fprintf(f,"\n");
        
            %%% harmonic bond
            elseif pot.style == "harmonic"
                fprintf(f,strcat(" ",...
                    num2str(pot.params.k_x/2)," ",...
                    num2str(pot.r12_eq),"\n"));
        
            %%% cut harmonic bond
            elseif pot.style == "harmonic_cut"
                fprintf(f,strcat(" ",...
                    num2str(pot.params.U_cut)," ",...
                    num2str(pot.r12_eq)," ",...
                    num2str(pot.params.r12_cut)," ","\n"));

            end
        end

    end
end