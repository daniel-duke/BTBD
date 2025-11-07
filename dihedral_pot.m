%%% dihedral potential class for BTBD
classdef dihedral_pot
    properties
        label           % name
        style           % dihedral style
        phi_eq          % equilibrium dihedral
        params          % parameters object
        index           % dihedral type (numeric)
    end

    methods
        %%% constructor
        function dpot = dihedral_pot(label,phi_eq,params_input,index)
            if nargin > 0
                dpot.label = label;
                dpot.phi_eq = phi_eq;
                dpot.index = index;
                dpot = dpot.set_params(params_input);
            end
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%% put names to parameters, convert energy to kcal/mol
        function dpot = set_params(dpot,params_input)

            %%% harmonic dihedral
            if dpot.label == "harmonic"
                dpot.style = "quadratic";
                dpot.params.k_phi = 6.96*params_input(1);
                if dpot.phi_eq > 90
                    warning("Equilibrium angle should be acute for quadratic dihedral to avoid instabilities.")
                end

            %%% cosine dihedral
            elseif dpot.label == "cosine"
                dpot.style = "harmonic";
                dpot.params.k_phi = 6.96*params_input(1);

                %%% set equilibrium angle parameter
                if dpot.phi_eq == 0
                    dpot.params.d = -1;
                elseif dpot.phi_eq == 180
                    dpot.params.d = 1;
                else
                    error("Equilibrium angle must be 0 or 180 for cosine dihedral.")
                end

            %%% error
            else
                error("Unknown dihedral type.")
            end
        end


        %%% calculate energy at separation distance
        function U = calc_energy(dpot,phi,phi_ref)
            arguments
                dpot; phi
                phi_ref = dpot.phi_eq
            end

            %%% harmonic dihedral
            if dpot.label == "harmonic"
                dphi = ars.applyPBC(phi-phi_ref,360);
                U = dpot.params.k_phi/2*dphi^2;

            %%% cosine dihedral
            elseif dpot.label == "cosine"
                dphi = ars.applyPBC(phi-phi_ref,360);
                U = dpot.params.k_phi*(1-cosd(dphi));

            end
        end
        
        
        %%% write coeffecients for potential
        function write_potential(dpot,f,is_hybrid)
            fprintf(f,strcat(...
                "dihedral_coeff  ",...
                num2str(dpot.index)));
            if is_hybrid
                fprintf(f," " + dpot.style);
            end

            %%% harmonic dihedral
            if dpot.label == "harmonic"
                fprintf(f,strcat(" ",...
                    num2str(dpot.params.k_phi/2)," ",...
                    num2str(dpot.phi_eq),"\n"));
            
            %%% cosine dihedral
            elseif dpot.label == "cosine"
                fprintf(f,strcat(" ",...
                    num2str(dpot.params.k_phi)," ",...
                    num2str(dpot.params.d)," 1\n"));

            end
        end

    end
end