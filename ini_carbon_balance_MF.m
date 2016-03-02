function [rlc_pt1,rrc_pt1,leaf_res,root_res,rgtot,rmtot,totalC,Assim] = ...
    ini_carbon_balance_MF(Tleaf,CO2,PAR,sunrise,sunset,rsratio,rosette_area,leaf_c,root_c,suc_c,sta_c,MF_c,timestep,growth_capacity,p)

%Initial values
%______________

vlmax25 = p(37); %Initial photosynthetic rubisco capacity per unit leaf area at 25degC (micromol CO2 m-2 s-1)


%Initial setting/conditions
%__________________________

daylength = sunset - sunrise;
        
is_light = 1;            
    
        
%Calculating photosynthesis
%__________________________
        
[net_rate] = photosynthesis(CO2,Tleaf,PAR,vlmax25,daylength,p);
        
        
%Calculating maintenance respiration
%___________________________________
        
[leaf_res,root_res] = mainres(Tleaf,leaf_c,root_c,suc_c,rosette_area,timestep,p);
        
        
%Calculating Assimilatory flux
%_____________________________

sta_c_endday = sta_c;
MF_c_endday = MF_c;

dummy_starch_consumption = 0;% doesn't get used anyway because is_light=1

[suc_sta_base,suc_MF_base,sta_use,MF_use,suc_equi,al_suc,suc_c_disp,suc_c_interm,Assim]...
= assimilation_MF_clock(daylength,is_light,net_rate,timestep,leaf_res,...
  root_res,suc_c,rosette_area,sta_c_endday,MF_c_endday,leaf_c,dummy_starch_consumption,p);


%Calculating organ demand
%________________________

[rrc_pt,totdem,rlc_pt,root_growth,leaf_growth] = organdemand(timestep,rsratio,leaf_c,growth_capacity,p);


%Calculating allocation 
%______________________

[rlc_pt1,suc_sta,root_gro1,rrc_pt1,leaf_gro1]...
= allocation(rrc_pt,totdem,rlc_pt,root_growth,leaf_growth,suc_c_disp,is_light);


%Calculating translocation
%_________________________

[leaf_trans,root_trans] = translocation(rosette_area,suc_c_interm,suc_equi,leaf_c,root_c);

%Calculating amount in each carbon pool
%______________________________________


rgtot = rlc_pt1 + rrc_pt1; %Total growth respiration: gC/plant
rmtot = leaf_res + root_res; %Total maintenance respiration: gC/plant 
        
totalC= leaf_c + root_c + rmtot + rgtot + sta_c + suc_c; %total C (g/plant)



