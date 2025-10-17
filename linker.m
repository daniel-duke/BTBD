%%% linker class for BTBD
classdef linker
    properties
        label           % linker name
        pot_index       % bond type
        r12_cut         % cutoff radius
        apot_index      % angle type
        theta_min       % minimum theta
        theta_max       % maximum theta
        ti_start        % starting atom type
        nlink5          % number of 5p sites
        link5s_oi       % 5p sites origami index
        link5s_bi       % 5p sites block index
        link5s_ib       % 5p sites patch bead index
        nlink3          % number of 3p sites
        link3s_oi       % 3p sites origami index
        link3s_bi       % 3p sites block index
        link3s_ib       % 3p sites patch bead index
    end

    methods
        %%% constructor
        function l = linker(label,pot_index,r12_cut,ti_start)
            if nargin > 0
                l.label = label;
                l.pot_index = pot_index;
                l.r12_cut = r12_cut;
                l.apot_index = 0;
                l.theta_min = 0;
                l.theta_max = 180;
                l.ti_start = ti_start;
                l.nlink5 = 0;
                l.link5s_oi = [];
                l.link5s_bi = [];
                l.link5s_ib = [];
                l.nlink3 = 0;
                l.link3s_oi = [];
                l.link3s_bi = [];
                l.link3s_ib = [];
            end
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%% add linker site
        function l = add_link(l,is_5p,o,ois,bi,patch)

            %%% 5p site
            if is_5p

                %%% get block indices
                if strcmp(bi,'A')
                    bis = 1:length(o.bs);
                elseif strcmp(bi,'B')
                    bis = 1:length(o.bs)-1;
                else
                    bis = str2double(bi);
                end

                %%% set indices
                for oi = ois
                    for bi = bis
                        l.nlink5 = l.nlink5 + 1;
                        l.link5s_oi(l.nlink5) = oi;
                        l.link5s_bi(l.nlink5) = bi;
                        l.link5s_ib(l.nlink5) = o.bs(bi).get_ib_patch(patch);
                    end
                end
            
            %%% 3p site
            else

                %%% get block indices
                if strcmp(bi,'A')
                    bis = 1:length(o.bs);
                elseif strcmp(bi,'B')
                    bis = 1:length(o.bs)-1;
                else
                    bis = str2double(bi);
                end

                %%% set indices
                for oi = ois
                    for bi = bis
                        l.nlink3 = l.nlink3 + 1;
                        l.link3s_oi(l.nlink3) = oi;
                        l.link3s_bi(l.nlink3) = bi;
                        l.link3s_ib(l.nlink3) = o.bs(bi).get_ib_patch(patch);
                    end
                end

            end
        end


        %%% write fix bond create command
        function write_bond_create(l,f,react_every)
            ti1 = l.ti_start;
            ti2 = l.ti_start+1;
            fprintf(f,strcat(...
                "fix             ", l.label, " stickies bond/create"));
            if l.apot_index ~= 0
                fprintf(f,"/angle");
            end
            fprintf(f,strcat(...
                " ", num2str(react_every), " ",...
                num2str(ti1), " ",num2str(ti2), " ",...
                ars.fstring(l.r12_cut,0,2), " ", num2str(l.pot_index), " "));
            if l.apot_index ~= 0 
                fprintf(f,strcat(...
                    "atype ", num2str(l.apot_index), " ",...
                    "aconstrain ", ars.fstring(l.theta_min,0), " ",ars.fstring(l.theta_max,0), " "));
            end
            fprintf(f,strcat(...
                "iparam 1 ", num2str(ti1), " ",...
                "jparam 1 ", num2str(ti2), "\n"));
        end

    end
end