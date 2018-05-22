function [rrc_pt,totdem,rlc_pt,root_growth,leaf_growth] = organdemand(timestep,rsratio,leaf_c,growth_capacity,p)
%% calculates partitioning of growth carbon to roots, leaves, and respiration
% Inputs:
%    timestep - 1 h by model-wide, hard-wired default
%    rsratio - Ratio of root to shoot growth rate (dimensionless)
%    leaf_c - Leaf carbon (gC/plant)
%    growth_capacity - maximum growth
%    p - Vector of parameters
%
% Outputs:
%    rrc_pt - Root growth respiration In g C/plant/timestep
%    totdem - sum of all growth demand, including respiration
%    rlc_pt - Leaf growth respiration In g C/plant/timestep (required)
%    root_growth - Root growth in gC/plant/timestep
%    leaf_growth - Leaf growth in gC/plant/timestep
%
% Functions called:
%    none



%%

%Leaf growth
%___________

maxgrowth = min(growth_capacity,p(63)); %value from Simile: 17*0.001*24 or limited by protein synthesis;
leaf_growth = maxgrowth*leaf_c*timestep; %Maximal leaf growth (gC/plant/timestep)

rc_coef = p(64); %Growth respiration to total growth allocation

numrlc = leaf_growth*rc_coef;
denomrlc = 1-rc_coef;

rlc_pt = numrlc/denomrlc; % Growth respiration In g C/plant/timestep (required)


%Root growth
%___________

root_growth = leaf_growth*rsratio; %Root growth (g C/plant/time step)

numrrc = root_growth*rc_coef;
denomrrc = 1-rc_coef;

rrc_pt = numrrc/denomrrc; % Growth respiration In g C/plant/timestep (required)


%Total growth demand
%___________________

totdem = leaf_growth + root_growth + rlc_pt + rrc_pt;
