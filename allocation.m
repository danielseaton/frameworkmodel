function [rlc_pt1,suc_sta,root_gro1,rrc_pt1,leaf_gro1] = allocation(rrc_pt,totdem,rlc_pt,root_growth,leaf_growth,suc_c_disp,is_light)
%% allocation of carbon to growth and respiration
% Input:
%   rrc_pt - Root growth respiration In g C/plant/timestep
%   totdem - sum of all growth demand, including respiration
%   rlc_pt - Leaf growth respiration In g C/plant/timestep (required)
%   root_growth - Root growth in gC/plant/timestep
%   leaf_growth - Leaf growth in gC/plant/timestep
%   suc_c_disp - Amount of sugar available for growth
%   is_light - binary, whether it is light or dark
%   
%
% Output:
%   suc_growth - sucrose for growth (g C/plant)
%   rrc_ptl - Root growth respiration In g C/plant/timestep (actually used)
%   root_gro1 - root growth actually achieved
%   rlc_ptl - Leaf growth respiration In g C/plant/timestep (actually used)
%   leaf_gro1 - leaf growth actually achieved

%%

if      totdem < suc_c_disp
    
        suc_growth = totdem;
        rrc_pt1 = rrc_pt;
        root_gro1 = root_growth;
        rlc_pt1 = rlc_pt;
        leaf_gro1 = leaf_growth;
        
        if  is_light == 1
            suc_sta = suc_c_disp - suc_growth;
        else
            suc_sta = 0;
        end    

else
        suc_growth = suc_c_disp; %sucrose for growth (g C/plant)
        rrc_pt1 = rrc_pt*(suc_growth/totdem); % Root growth respiration In g C/plant/timestep (actually used)
        root_gro1 = root_growth*(suc_growth/totdem); %root growth actually achieved
        rlc_pt1 = rlc_pt*(suc_growth/totdem); % Leaf growth respiration In g C/plant/timestep (actually used)
        leaf_gro1 = leaf_growth*(suc_growth/totdem); %leaf growth actually achieved (g C/plant/timestep)
        suc_sta = 0;
end


