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

            %%% cosine harmonic dihedral
            if dpot.label == "harmonic_cos"
                dpot.style = "harmonic";
                dpot.params.k_phi = 6.96*params_input(1);

                %%% set equilibrium angle parameter
                if dpot.phi_eq == 0
                    dpot.params.d = 1;
                elseif dpot.phi_eq == 180
                    dpot.params.d = -1;
                else
                    error("Equilibrium angle must be 0 or 180 for harmonic cosine dihedral.")
                end

            %%% quadratic harmonic dihedral
            elseif dpot.label == "harmonic_quad"
                dpot.style = "quadratic";
                dpot.params.k_phi = 6.96*params_input(1);

            %%% error
            else
                error("Unknown dihedral type.")
            end
        end


        %%% calculate energy at separation distance
        function U = calc_energy(dpot,phi)

            %%% cosine harmonic dihedral
            if dpot.label == "harmonic_cos"
                U = dpot.params.k_phi*(1-dpot.params.d*cosd(phi));

            %%% quadratic harmonic dihedral
            elseif dpot.label == "harmonic_quad"
                U = dpot.params.k_phi/2*(phi-dpot.params.phi_eq)^2;

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

            %%% cosine harmonic dihedral
            if dpot.label == "harmonic_cos"
                fprintf(f,strcat(" ",...
                    num2str(dpot.params.k_phi/2)," ",...
                    num2str(dpot.params.d)," 1\n"));

            %%% cosine harmonic dihedral
            elseif dpot.label == "harmonic_quad"
                fprintf(f,strcat(" ",...
                    num2str(dpot.params.k_phi)," ",...
                    num2str(dpot.phi_eq),"\n"));
            
            end
        end

    end
end