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
            end
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%% calculate separation distance at energy
        function r12 = calc_separation(pot,U)
        
            %%% harmonic bond
            if pot.style == "harmonic"
                k_x = pot.params(1);
                r12 = pot.r12_eq + sqrt(2*U/k_x);
        
            %%% cut harmonic bond
            elseif pot.style == "harmonic/shift/cut"
                k_x = pot.params(1);
                U_cut = pot.params(2);
                r12_Ucut = pot.r12_eq + sqrt(2*U_cut/k_x);
                r12_U = pot.r12_eq + sqrt(2*U/k_x);
                r12 = min([r12_Ucut,r12_U]);
        
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
                pot.style, " "));
        
            %%% harmonic bond
            if pot.style == "harmonic"
                k_x = pot.params(1);
                fprintf(f,strcat(...
                    num2str(k_x)," ",...
                    num2str(pot.r12_eq),"\n"));
        
            %%% cut harmonic bond
            elseif pot.style == "harmonic/shift/cut"
                k_x = pot.params(1);
                U_cut = pot.params(2);
                r12_cut = pot.r12_eq + sqrt(2*U_cut/k_x);
                fprintf(f,strcat(...
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