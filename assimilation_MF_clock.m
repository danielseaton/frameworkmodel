function [suc_sta_base,suc_MF_base,sta_use,MF_use,suc_equi,al_suc,suc_c_disp,suc_c_interm,al_pt_plant_assim]...
    = assimilation_MF_clock(daylength,is_light,net_rate,timestep,leaf_res,...
    root_res,suc_c_perplant,rosette_area,sta_c_endday,MF_c_endday,leaf_c,starch_consumption,p,mf_use)

convert_to_gC = timestep*p(59)*10^(-6)*12; %conversion factor for umol/m2 leaf/sec to gC/m2 leaf/timestep

%Baseline conversion coefficient (default p(60)=0.125)
%sta_base = p(60); %Baseline starch conversion coefficient

%Baseline conversion coefficient based on TiMet's AGPase activity data
sta_base = 0.6*(-0.0296*daylength + 0.7157);

%Malate and fumarate conversion coefficient based on Laurel and Hardy Jan
MF_base = 0.2*0.6*(-0.0296*daylength + 0.7157);
%MF_base = sta_base*(-0.0051*daylength + 0.1278); %for WP12A
  
%sta_convert_night = p(61); % DANIEL: this parameter is now unused
%sta_convert_night = 0.89;
MF_convert_night = mf_use;
%MF_convert_night = -0.0325*daylength + 1.1033; %for WP12A

if      is_light == 1
    
        sta_use = 0; %conversion of starch to sugar
        MF_use = 0; %conversion of fumarate and malate to sugar
        al_pt_plant_assim = net_rate*rosette_area*convert_to_gC; %Assimilatory flux per plant (gC/plant/timestep)
        suc_sta_base = al_pt_plant_assim*sta_base; %baseline starch conversion
        suc_MF_base = al_pt_plant_assim*MF_base; %baseline fm conversion
        al_suc =al_pt_plant_assim - suc_sta_base - suc_MF_base;
else
        al_pt_plant_assim = 0;
        sta_use = starch_consumption; %Dark conversion of starch to sugar. DANIEL: really simple change here
        %sta_use = ((sta_c_endday*sta_convert_night)/(24-daylength))*timestep*24;
        MF_use = ((MF_c_endday*MF_convert_night)/(24-daylength))*timestep*24;
        suc_sta_base = 0;
        suc_MF_base = 0;
        al_suc = 0;
end

suc_equi = p(62); %equilibrium sucrose plus hexose concentration in leaves (g C/m2 leaf)
%suc_equi = 0.000315; %lowest measured sucrose level in pgm from Pal et al 2013 (0.17 micromol suc/g FW = 0.000315 gC/g DW based on Col water content of 92%)

rosette_DW = leaf_c/p(30); %rosette DW calculated from leaf C and carbon content
suc_c_interm = suc_c_perplant + sta_use + MF_use + al_suc - root_res - leaf_res;

%Growth capacity limited by polysome loading (a proxy of protein synthesis rate) 
%growth_capacity = 92.681*suc_c_perplant/rosette_DW + 0.0958; %linear regression from Pal et al 2013 for Col wt

%suc_equi = (sta_use + al_suc - root_res - leaf_res)/rosette_DW*(0.5 - growth_capacity/24);
%suc_equi = max(0,suc_equi);
%suc_equi = (0.1288*daylength+0.1157)*12*1e-06*12/(1-0.92); %minimum suc at dawn from WP1.2A

%current_value = suc_c_interm - (suc_equi*rosette_DW);
current_value = suc_c_interm - (suc_equi*rosette_area);


if      current_value <= 0
        suc_c_disp = 0;
else
        suc_c_disp = current_value; %Amount of sugar available for growth
end

%growth_capacity = 1;
                                                            
                                                         
