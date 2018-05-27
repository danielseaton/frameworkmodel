function [leaf_trans,root_trans] = translocation(rosette_area,suc_c_interm,suc_equi,leaf_c,root_c)
%% transolcation of sugar between compartments if sucrose falls below equilibrium values
%
% Input:
%   rosette_area - Total leaf area of the rosette
%   suc_c_interm - Intermediate sucrose value
%   suc_equi - equilibrium sucrose plus hexose concentration in leaves (g C/m2 leaf)
%   leaf_c - Leaf carbon (gC/plant)
%   root_c - Root carbon (gC/plant)
%
% Output:
%   leaf_trans - Translocation from leaves (gC/plant)
%   root_trans - Translocation from roots (gC/plant)


global p

suc_equi_plant = suc_equi*rosette_area; %Equilibrium sucrose plus hexose mass for whole plant (g C/plant)
%suc_equi_plant = suc_equi*leaf_c/p(30); %Equilibrium sucrose plus hexose mass for whole plant (g C/plant)
root_and_leaf_c = root_c + leaf_c; %Total root and leaf C (g C/plant)

if      suc_c_interm <= suc_equi_plant %translocation needed
        
        scalingl = leaf_c/root_and_leaf_c;
        leaf_trans = (suc_equi_plant - suc_c_interm)*scalingl; %Translocation from leaves (gC/plant)
        scalingr = root_c/root_and_leaf_c;
        root_trans = (suc_equi_plant - suc_c_interm)*scalingr; %Translocation from root (g C/plant)
        
else
        leaf_trans = 0;
        root_trans = 0;
end        
    
