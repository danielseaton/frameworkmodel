function [rlc_pt1,rrc_pt1,leaf_res,root_res,leaf_carbon,root_carbon,...
         sucrose_carbon,starch_carbon,MF_carbon,rgtotal,rmtotal,totalCarbon,Assim,suc_sta,net_rate,new_starch_module_state]...
         = plant_carbon_balance_MF_clock(hour_idx,Tleaf,CO2,PAR,sunrise,sunset,is_light,...
           rsratio,rosette_area,leaf_c,root_c,suc_c,sta_c,MF_c,rgtot,rmtot,...
           timestep,sta_c_endday,MF_c_endday,efficiency,growth_capacity,...
           clock_output,starch_module_state,starch_parameters,p,mf_use)

%Initial values
%______________

vlmax20 = p(38); %Photosynthetic rubisco capacity per unit leaf area at 20degC (micromol CO2 m-2 s-1)
vlmax25 = vlmax20/0.64;

%Initial setting
%_______________

daylength = sunset - sunrise;
      
     
        
%Calculating photosynthesis
%__________________________
        
if      is_light == 1
              
            [net_rate] = efficiency*photosynthesis(CO2,Tleaf,PAR,vlmax25,daylength,p); %micromol CO2 m-2 leaf s-1        
              
else
        net_rate = 0;
end    
        
        
%Calculating maintenance respiration
%___________________________________
        
[leaf_res,root_res] = mainres(Tleaf,leaf_c,root_c,suc_c,rosette_area,timestep,p);

% DANIEL
%Calculating starch degradation. Start by just connecting it in a trivial
%way, reading out starch content and regulating X appropriately.
[new_starch_module_state,starch_consumption] = starch_module(hour_idx,starch_module_state,sta_c,leaf_c,is_light,clock_output,starch_parameters);

%[starch_consumption] =
%starch_module(t,clock_output,sta_c,leaf_c,is_light);




%Calculating Assimilatory flux
%_____________________________     

%net_rate is the rate of photosynthesis
[suc_sta_base,suc_MF_base,sta_use,MF_use,suc_equi,al_suc,suc_c_disp,suc_c_interm,Assim]...
= assimilation_MF_clock(daylength,is_light,net_rate,timestep,leaf_res,...
root_res,suc_c,rosette_area,sta_c_endday,MF_c_endday,leaf_c,starch_consumption,p,mf_use);
      
      
%Calculating organ demand
%________________________    
	
[rrc_pt,totdem,rlc_pt,root_growth,leaf_growth] = organdemand(timestep,rsratio,leaf_c,growth_capacity,p);



%Calculating allocation 
%______________________

[rlc_pt1,suc_sta,root_gro1,rrc_pt1,leaf_gro1]...
= allocation(rrc_pt,totdem,rlc_pt,root_growth,leaf_growth,suc_c_disp,is_light);
%suc_sta which represented the 'overflow' to starch is no longer used

%Calculating translocation
%_________________________

[leaf_trans,root_trans] = translocation(rosette_area,suc_c_interm,suc_equi,leaf_c,root_c);

        
%Calculating amount in each carbon pool
%______________________________________

leaf_carbon = leaf_c + leaf_gro1 - leaf_trans; % Leaf carbon: g C/plant
starch_carbon = sta_c + suc_sta_base - sta_use; %Starch content per plant: gC/plant
MF_carbon = MF_c + suc_MF_base - MF_use; %Starch content per plant: gC/plant
root_carbon = root_c + root_gro1 - root_trans; %Root carbon: gC/plant
sucrose_carbon = suc_c - rlc_pt1 - rrc_pt1 - leaf_res...
                 - root_res - root_gro1 + root_trans...
                 + al_suc - leaf_gro1 + leaf_trans...
                 + sta_use + MF_use; %Sucrose carbon: gC/plant
rgtotal = rgtot + rlc_pt1 + rrc_pt1; %Total growth respiration: gC/plant
rmtotal = rmtot + leaf_res + root_res; %Total maintenance respiration: gC/plant 
        
totalCarbon= leaf_carbon + root_carbon + rmtotal + rgtotal + starch_carbon + MF_carbon + sucrose_carbon; %total C (g/plant)


        

