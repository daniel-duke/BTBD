%%% bond potential class for BTBD
classdef reaction
    properties
        style           % reaction style
        params          % other parameters
        index           % reaction number
    end

    methods
        %%% constructor
        function rxn = reaction(style,params,index)
            if nargin > 0
                rxn.style = style;
                rxn.params = params;
                rxn.index = index;
            end
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        

    end
end