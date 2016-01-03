% function [cvRMSE] = integrationrebalMF_clock(temp,rise,set,co2,light)
%[cvRMSE]=integrationrebalMF_clock(20.5,0,12,42,145)


clock_genotype = {''};
% clock_genotype = {'prr9','prr7'};


temp = 20.5;
rise = 0;
set = 12;
co2 = 42;
light = 145;

%specify Col accession for phenology model threshold
geno = 2;

% %Mugford
% temp = 20;
% rise = 0;
% set = 12;
% co2 = 42;
% light = 190;

% %Ishihara
% temp = 19.3333;
% rise = 0;
% set = 8;
% co2 = 42;
% light = 150;

% %Comparot-Moss conditions
% temp = 20;
% rise = 0;
% set = 12;
% co2 = 42;
% light = 170;

genotype = 1;

%Calling for parameters
%______________________

global p d mf_use

addpath('PIF_CO_FT_model')

load('parameter.mat')

p=parameter;

probability = 1; %deterministic

switch genotype
    case 1
        w = 0.91; %water content
        mf_use = 0.7; %malate+fumarate turnover
    case 2
        w = 0.89; %water content
        mf_use = 0.7; %malate+fumarate turnover
    case 3
        w = 0.89; %water content
        mf_use = 0.7; %malate+fumarate turnover
    case 4
        w = 0.89; %water content
        mf_use = 0.25; %malate+fumarate turnover
end        

% w = 0.9;

d = 1-w; %dry matter


%Calling for meteorological data
%_______________________________  
    
    

load('weather.mat')

hour=weather(:,1);
T=temp*weather(:,2);
sunrise=rise*weather(:,3);
sunset=set*weather(:,4); 
CO2=co2*weather(:,5); %CO2 partial pressure (Pa)
PAR=light*weather(:,6); %total absorbed PAR per unit leaf area (micromol m-2 s-1)
Photoperiod = sunset(end) - sunrise(end);


%Parameters for growth model
%__________________________________________

Tb = p(1); %Base temperature for thermal-time calculation
TT0 = p(19); %degree days required from seed germination until plant emergence (stage 1.0 in Boyes et al. 2001)


TS = p(26); %degree days of leaf lifespan (from the point of leaf initiation)
TLcot = p(27); %cotyledon expansion period
TLtrue = p(28); %expansion period of true leaves
Root_LC = p(29); %root expansion period 
leaf_factor = p(30); %carbon-dry mass conversion
root_factor = p(31); 
seed_input = p(32);
cot_area = p(33); %when they are fully opened
sucrose_initial = p(34); %initial sucrose content per unit area
starch_initial = p(35); %initial starch content per initial sucrose content
petiole = p(36); %maximum petiole elongation factor 
min_alpha = p(65); %minimum zenithal angle
inc_alpha = p(66); %maximum increase in zenithal angle from the minimum value
juv_TT = p(67); %thermal time since plant emergence for the juvenile-to-adult transition
juv_rate = p(68); %rate of leaf appearance at juvenile stage
juv_interc = p(81);%y-intercept for the juvenile rate of leaf appearance
ad_rate = p(69); %rate of leaf appearance at adult stage
ad_interc = p(70); %y-intercept for the adult rate of leaf appearance
SLA_cot = p(71); %SLA for cotyledons (in m2/g dry mass) Fig. 3 Christophe et al (2008)
SLA_exp = p(72); %curve shape for SLA against thermal time

clear starch
clear sugar
clear GPP
clear NPP
clear GPP_day
clear NPP_day 
clear net_to_gross 
clear net_to_gross_day 
clear growth_res_ratio 
clear Gas_exchange_perFW 
clear starch_perFW 
clear sugar_perFW 



%%%This is parameterisation - parameters derived from other parameters but
%%%not depending on conditions
%Expansion period for the root system
%____________________________________

PR=p(23); %root expansion
aR=p(24);
bR=p(25);

TR = 1000;%Root_LC*CumThrm(end);
MaxfRpoint = TR/(1+(bR-1)/(aR-1))-0.5; %determining the maximum point by solving the function derivative
MaxfR = ((MaxfRpoint+0.5)/TR)^(aR-1)*(1-(MaxfRpoint+0.5)/TR)^(bR-1);  

%Maximum values of leaf sink variation
%_____________________________________

PL=p(20); %leaf expansion
aL=p(21);
bL=p(22);

Maxfcotpoint = TLcot/(1+(bL-1)/(aL-1))-0.5; %determining the maximum point by solving the function derivative
Maxfcot = ((Maxfcotpoint+0.5)/TLcot)^(aL-1)*(1-(Maxfcotpoint+0.5)/TLcot)^(bL-1);

MaxfLpoint = TLtrue/(1+(bL-1)/(aL-1))-0.5; %determining the maximum point by solving the function derivative
MaxfL = ((MaxfLpoint+0.5)/TLtrue)^(aL-1)*(1-(MaxfLpoint+0.5)/TLtrue)^(bL-1);  
          
%Initialise matrices that will fill up
fL = [];
fR = [];
rsratio = [];
efficiency = [];
growth_capacity = [];



%%%This should be a separate function
%initialise clock to starting conditions
clock_state_0=[1.0151 0.956 0.0755 0.0041 0.506 0.0977 0.0238 0.0731 0.0697 0.0196 0.0435 0.2505 0.0709 0.1017 0.0658 0.4016 0.1167 0.1012 0.207 0.0788 0.3102 0.0553 0.2991 0.1503 0.0286 0.65 0.2566 0.1012 0.576 0.3269]; %12:12 wt;
for i = 1:5
    clock_output=circadian_module(sunrise(1),sunset(1),clock_state_0,clock_genotype);
    %work out the clock state at ZT24 i.e. at the end of the previous day
    clock_state_0 = interp1q(clock_output.T,clock_output.Y,24);
end

starch_module_state = [1,1,1,1];
FT_module_state =ones(1,18);

% clock_genotype = {''};

% set whether the phenology model controls the end of the simulation
N_max_days = 40; %days max simulation
input_size = 40*24;
run_phenology_model = 0;

day_idx = 1;
has_flowered = false;
has_emerged = false;
CumPhenThrm=0;

while day_idx <= N_max_days && ~(has_flowered)
    %initial timepoint
    t = (day_idx-1)*24+1;
    
    %run clock model for this day
    clock_output=circadian_module(sunrise(t),sunset(t),clock_state_0,clock_genotype);
    %work out the clock state at ZT24 i.e. at the end of the previous day
    clock_state_0 = interp1q(clock_output.T,clock_output.Y,24);
    
    %run phenology model
    if run_phenology_model
        [DayPhenThrm,FT_module_state] = phen(T,sunrise,sunset,geno,clock_output,FT_module_state);
    end
    DayPhenThrm = 0;
    CumPhenThrm = DayPhenThrm+CumPhenThrm;
    has_flowered = flowering_threshold_test(CumPhenThrm,genotype);
    
    for hour_idx = 1:24
        %timepoint
        t = (day_idx-1)*24+hour_idx;
        
        Thrm(t)=max(0,(T(t)-Tb))/24; %in degree days   
        CumThrm(t) = sum(Thrm(1:t)); %cumulative thermal time
        
        if ~(has_emerged)
            %Plant has not yet emerged
            %______________________________________________

            if CumThrm(t) >= TT0
                eme = t;
                has_emerged = 1;
                %For the first growth cycle (from germination to TT0):
                %____________________________________________________
                %%%FUNC: initialise seedling state at emergence

                Seed_mass = seed_input; %Seed dry mass (g), emptied after the first cycle to reach stage 1.0
                Q(eme) = Seed_mass;
                fR = ((CumThrm(eme)+0.5)/TR)^(aR-1)*(1-(CumThrm(eme)+0.5)/TR)^(bR-1)/MaxfR;  
                D(eme) = PR*fR + 2*PL;
                Root_mass(eme) = PR*fR*Q(eme)/D(eme); %root dry mass (g)

                Si(eme,1:2) = cot_area; %area for EACH cotyledon (m2), cycle 1 for leaf ranked 1 and 2
                Leaf_mass(eme,1:2) = Si(eme,1)/SLA_cot; %dry mass for EACH cotyledon (g)
                Hypocotyl_mass = Seed_mass - Root_mass(eme) - sum(Leaf_mass(eme,:));
                Total_shoot(eme) = Seed_mass - Root_mass(eme);
                Destructive_area(eme) = sum(Si(eme,:));
                Leaf_no(eme) = 2;
                Appear(1:2) = eme; %The timepoint when the cotyledons start expanding using photosynthate


                %Initialisation
                %______________

                S(eme) = sum(Si(eme,1:2));
                S_intercept(eme) = S(eme)*cos(10/180*pi); 
                Func_area(eme)=S(eme);
                %%%This could give the impression that another timestep can
                %%%be chosen
                steps_perhour = 1;
                steps_perday = 24*steps_perhour;
                timestep = 1/steps_perday; % Time step length used to solve the model

                Leaf_carbon(eme) = leaf_factor*sum(Leaf_mass(eme,:));
                Root_carbon(eme) = root_factor*Root_mass(eme);    
                Sucrose_carbon(eme) = sucrose_initial*S(eme); %Initial amount of sucrose carbon (g/plant)
                Starch_carbon(eme) = starch_initial*Sucrose_carbon(eme); %Initial amount of starch carbon (g/plant)
                MF_carbon(eme) = 0.4*Starch_carbon(eme); %Initial amount of malate+fumarate carbon (g/plant)
                rsratio(eme) = Root_carbon(eme)/Leaf_carbon(eme);
                growth_capacity(eme) = p(63);

                if      hour(eme) < sunrise(eme)
                        sta_c_endday = Starch_carbon(eme)/(p(61)*((sunrise(eme)-hour(eme))/(sunrise(eme)+24-sunset(eme-1)))+1-p(61));
                elseif  hour(eme) > sunset(eme)
                        sta_c_endday = Starch_carbon(eme)/(p(61)*((24-hour(eme)+sunrise(eme+1))/(24-sunset(eme)+sunrise(eme+1)))+1-p(61));
                else    
                        sta_c_endday = 0;
                end        

                MF_c_endday = 0.4*sta_c_endday;

                [rlc_pt1(eme),rrc_pt1(eme),leaf_res(eme),root_res(eme),rgtot(eme),rmtot(eme),totalC(eme),assim(eme)] = ...
                                                      ini_carbon_balance_MF(T(eme),CO2(eme),PAR(eme),...
                                                      sunrise(eme),sunset(eme),rsratio(eme),...
                                                      S_intercept(eme),Leaf_carbon(eme),...
                                                      Root_carbon(eme),Sucrose_carbon(eme),...
                                                      Starch_carbon(eme),MF_carbon(eme),timestep,growth_capacity(eme));
                %%%END FUNC: initialise seedling state at emergence

                GC = 1; %number of growth cycle elapsed                                  
                GCstart(GC+1) = eme; %The timepoint when the second growth cycle starts
            end
        else
            %Plant is growing
            %________________

            %To determine the rosette area at the previous time point for light interception
            %_________________________________________________________________________________


            %Determining for each leaf rank i (functioning/senesced) at the previous timepoint:

            for i = 1:Leaf_no(t-1)
                if CumThrm(t-1) < (CumThrm(Appear(i))+TS)
                   fS(i) = 1;
                else
                   fS(i) = 0;
                end
            end


            %To determine rosette structure (at the previous time point):

            [Max_area,i_max] = max(Si(t-1,:));
            i_max = max(i_max,2);

            %%%DS note - change loop or remove mult by indicator variable fS(i)
            for i=1:Leaf_no(t-1) %for each leaf rank i

                if      i<=i_max

                        alpha(i)=min_alpha; %zenithal angle
                else 
                        alpha(i)=min_alpha+inc_alpha*(i-i_max)/(Leaf_no(t-1)-i_max); 
                end
               
                S(t-1,i)=Si(t-1,i)*cos(alpha(i)/180*pi)*fS(i); %Projected functioning leaf area 
                Sf(t-1,i)=Si(t-1,i)*fS(i); %Functioning leaf area 

            end


            %Rosette area (at the previous time point):      
            %FUNC: calculate area a previous time point
            [S_intercept,Func_area] = calculate_area(t,S,petiole,input_size,Appear,Leaf_no,Sf,S_intercept,Func_area);
            %END FUNC:  calculate area a previous time point


            %To determine whole-plant carbon balance (at the current time point):
            %___________________________________________________________________


            if CumThrm(t) < (juv_TT+TT0)
               phyllochron = 30.3*probability;
            else
               phyllochron = 11.9*probability;
            end

            %To determine current leaf number:
            %FUNC: determine current leaf number
            [Appear,Leaf_mass,Si,Leaf_no,GC,GCstart] = determine_current_leaf_number(...
                t,probability,CumThrm,phyllochron,Appear,Leaf_mass,Si,Leaf_no,GC,GCstart);
            %END FUNC: determine current leaf number

            %FUNC: root shoot partitioning
            %To calculate sink variation for leaf:
            [fL,fR,rsratio] = calculate_root_shoot_sink_variation(t,Leaf_no,TLcot,Maxfcot,...
                TLtrue,MaxfL,aL,bL,CumThrm,Appear,TR,MaxfR,aR,bR,PR,PL,rsratio,...
                leaf_factor,root_factor,fL,fR);
            %Calculating root-to-shoot ratio (in g Carbon):
            leafconversion = leaf_factor;     % Gorsuch et al 2010a and 2010b

            rootconversion = root_factor; %g Carbon per g dry mass
                             %from Kumar et al (Agroforest Syst 2010) around 30-35%
                             %from UK BioChar Research Centre (UoE) around 35 % in
                             %wheat, cited as Prendergast-Miller M and Sohi SP 2010. 
                             %Investigating biochar impacts on plant roots and root carbon.
                             %Poster presented in the Organic Matter Stabilization
                             %and Ecosystem Functions session at Soil Organic
                             %Matter Conference, Cote d'Azur, France (Sept 2010)

            % To determine light status:                
            %islight = determine_light_status(t,hour,sunrise,sunset)
            if  hour(t) > sunrise(t) && hour(t) <= sunset(t)            
                is_light(t) = 1;            
            else
                is_light(t) = 0;
            end   


            %To determine starch status:
            %sta_c_endday = determine_starch_status(t,sta_c_endday,Starch_carbon,is_light)
            if      is_light(t) == 1     
                    sta_c_endday = 0;

            elseif  is_light(t) == 0 && is_light(t-1) == 1  
                    sta_c_endday = Starch_carbon(t-1);
            else   
                    sta_c_endday = sta_c_endday;
            end


            %To determine malate and fumarate (MF) status:
            %MF_c_endday = determine_MF_status(t,MF_c_endday,MF_carbon,is_light)

            if      is_light(t) == 1     
                    MF_c_endday = 0;

            elseif  is_light(t) == 0 && is_light(t-1) == 1  
                    MF_c_endday = MF_carbon(t-1);
            else   
                    MF_c_endday = MF_c_endday;
            end


            %To determine the amount of sugar accumulated over a day
            %%FUNC: [efficiency,growth_capacity]=calculate_sugar_acc(t,eme,efficiency,growth_capacity,is_light)
            [efficiency,growth_capacity]=calculate_sugar_acc(t,p,eme,...
                efficiency,growth_capacity,is_light);
            %END FUNC: calculate_sugar_acc

            %DANIEL: added new_clock_state output and clock_state input
            [rlc_pt1(t),rrc_pt1(t),leaf_res(t),root_res(t),...
             Leaf_carbon(t),Root_carbon(t),Sucrose_carbon(t),Starch_carbon(t),MF_carbon(t),rgtot(t),rmtot(t),totalC(t),assim(t),...
             suc_sta(t),net_rate(t),new_starch_module_state] = plant_carbon_balance_MF_clock(hour_idx,T(t),CO2(t),PAR(t),sunrise(t),sunset(t),is_light(t),...
                   rsratio(t),S_intercept(t-1),Leaf_carbon(t-1),Root_carbon(t-1),...
                   Sucrose_carbon(t-1),Starch_carbon(t-1),MF_carbon(t-1),rgtot(t-1),rmtot(t-1),timestep,sta_c_endday,MF_c_endday,efficiency(t),...
                   growth_capacity(t),clock_output,starch_module_state);

            starch_module_state = new_starch_module_state;

            %Root mass:
            Root_mass(t) = Root_carbon(t)/rootconversion; %root dry mass (g)
            %Actual root C growth
            CarbonR = Root_carbon(t) - Root_carbon(t-1);    

            %Total dry biomass accumulation for leaves (g): 
            CarbonL = Leaf_carbon(t) - Leaf_carbon(t-1);
            BiomassL = CarbonL/leafconversion; 


            %Leaf demand:
            for     i = 1:Leaf_no(t)
                    Ld(t,i) = PL*fL(t,i);
            end           

            Totalleafdemand = sum(Ld(t,:));


            %To calculate SLA:        
            SLA(t) = SLA_cot*exp(-SLA_exp*(CumThrm(t)-TT0)); %SLA (in m2/g dry mass) Fig. 3 Christophe et al (2008)


            %Leaf mass and area:
            for     i = 1:Leaf_no(t)
                    Leaf_mass(t,i) = Leaf_mass(t-1,i) + Ld(t,i)/Totalleafdemand*BiomassL;
                    Si(t,i) = max(Leaf_mass(t,i)*SLA(t),Si(t-1,i));
            end        

            Destructive_area(t) = sum(Si(t,1:Leaf_no(t)));        
            Total_leaf_mass(t) = sum(Leaf_mass(t,:));
            Total_shoot(t) = Total_leaf_mass(t) + Hypocotyl_mass; 


            Cstatus(t) = (CarbonL+CarbonR)/(Totalleafdemand+fR(t)*PR);
            GroRes(t) = rgtot(t)-rgtot(t-1);
            MainRes(t) = rmtot(t)-rmtot(t-1);
            Total_plant(t) = Total_shoot(t) + Root_mass(t);  

            Normalisedgrowth(t) = (Total_plant(t)-Total_plant(t-1))/Total_plant(t-1);
            Normalisedcost(t) = net_rate(t)/Total_plant(t-1);

        end
    end
    
    
    %Update day index
    day_idx = day_idx+1;
end

%-------------------------------End-of-model-------------------------------
totgrowth = sum(Normalisedgrowth(eme+25:end));
totcost = sum(Normalisedcost(eme+25:end));
minratio = totcost/totgrowth;
Leaf_no(end)
S_intercept(end)
Root_mass(end)

Measured_day = [22 28 29 38];
Measured_hour = [0 12 24];

switch genotype
    case 1
        Measured_biomass = [0.034 0.1262 0.2192 0.95014]; %Col
        Measured_sd_biomass = [0.004092676	0.017887426	0.032882518	0.14407886]; %Col
        Measured_starch = [2.41 42.70 4.55]; %Col
        Measured_OA = [4.85 17.87 5.38]; %Col
        sd_starch = [0.35 2.37 0.31]; %Col
        sd_OA = [0.51 1.24 0.71]; %Col
        
    case 2
        Measured_biomass = [0.01866	0.08546	0.0968	0.48772]; %lsf
        Measured_sd_biomass = [0.003036939	0.007267255	0.012679511	0.046257616]; %lsf
        Measured_starch = [45.28 90.19 44.25]; %lsf
        Measured_OA = [6.62 24.18 7.85]; %lsf    
        sd_starch = [1.05 4.61 2.08]; %lsf
        sd_OA = [1.18 1.07 0.27]; %lsf
        
    case 3
        Measured_biomass = [0.02102	0.09238	0.10086	0.5743]; %prr7/prr9
        Measured_sd_biomass = [0.001413153	0.013399515	0.008260327	0.061144337]; %prr7/prr9
        Measured_starch = [9.93 67.31 23.77]; %prr7/prr9
        Measured_OA = [18.65 41.80 25.31]; %prr7/prr9    
        sd_starch = [0.08 0.01 0.14]; %prr7/prr9
        sd_OA = [2.14 2.69 4.04]; %prr7/prr9

    case 4
        Measured_biomass = [0.02102	0.09238	0.10086	0.5743]; %prr7/prr9
        Measured_sd_biomass = [0.001413153	0.013399515	0.008260327	0.061144337]; %prr7/prr9
        Measured_starch = [9.93 67.31 23.77]; %prr7/prr9
        Measured_OA = [18.65 41.80 25.31]; %prr7/prr9    
        sd_starch = [0.08 0.01 0.14]; %prr7/prr9
        sd_OA = [2.14 2.69 4.04]; %prr7/prr9
end        

Estimate_biomass = Total_shoot(Measured_day*24)/d;
cvRMSE = (sum((Estimate_biomass - Measured_biomass).^2)/numel(Measured_biomass))^0.5/...
                (mean(Measured_biomass));

nRMSE = (sum((Estimate_biomass - Measured_biomass).^2)/numel(Measured_biomass))^0.5/...
                (max(Measured_biomass) - min(Measured_biomass))
            
Rsquare = corrcoef(Measured_biomass,Estimate_biomass).^2;
Rsquare = Rsquare(2)
         
%sugar = [Sucrose_carbon(end)/Total_shoot(end) Sucrose_carbon(end-Photoperiod)/Total_shoot(end-Photoperiod)]*d/12/12*10^6; 

%contribution of starch and malate+fumarate
day = 1;
for x=252:24:input_size-12
    starchcontribution = Starch_carbon(x) - Starch_carbon(x+12);
    malatefumaratecontribution = MF_carbon(x) - MF_carbon(x+12);
    MFcontribution(day) = malatefumaratecontribution/(malatefumaratecontribution + starchcontribution)*100;
    day = day+1;
end

% figure
% plot(10:39,MFcontribution)
% axis([10 40 10 20])
% title('Percentage contribution of carbon from malate+fumarate at night')
% xlabel('DAS')
% ylabel('Percentage (% of total carbon degraded)')

op = 0;

if op == 0
    
Fresh_weight = Total_shoot(end)/d %g FW rosette
starch = Starch_carbon; %g C
sucrose = Sucrose_carbon; %g C
NPP = assim - rlc_pt1 - rrc_pt1 - leaf_res - root_res; %g C/h
Gas_exchange_perFW = NPP./(Total_shoot/d); %g C/g FW/h
starch_perFW = starch./(Total_shoot/d); %g C/g FW
sucrose_perFW = sucrose./(Total_shoot/d); %g C/g FW


FW_ED = Total_shoot(input_size+Photoperiod-24)/d; %g
FW_EN = Total_shoot(input_size)/d; %g
AD = Gas_exchange_perFW(input_size+2-24); %g C/g FW/h
AR1 = Gas_exchange_perFW(input_size-24+Photoperiod+5); %g C/g FW/h
AR2 = Gas_exchange_perFW(input_size-1); %g C/g FW/h


FW = Total_shoot/d;
starchC6 = Starch_carbon/12/6*10^6./FW;
malatefumarateC4 = MF_carbon/12/4*10^6./FW;
sugar = Sucrose_carbon/12/12*10^6./FW;
gFW_22 = FW(22*24)
gFW_28 = FW(28*24)
gFW_29 = FW(29*24)
gFW_38 = FW(38*24)
AperFW = Gas_exchange_perFW(38*24+1)/12*10^6 %micromol CO2/g FW/h
RperFW = Gas_exchange_perFW(38*24-1)/12*10^6 %micromol CO2/g FW/h
Aperarea = NPP(38*24+1)/S_intercept(38*24)/12*10^6/10000 %micromol CO2/cm2/h
Rperarea = NPP(38*24-1)/S_intercept(38*24-2)/12*10^6/10000 %micromol CO2/cm2/h


figure
hold on
plot(0:24,starchC6(28*24:29*24))
errorbar(Measured_hour,Measured_starch,sd_starch,'o')
axis([0 24 0 100])
title('Starch profile')
xlabel('Hour')
ylabel('Starch level (micromol C6 per g FW)')
hold off

% figure
% hold on
% plot(0:24,malatefumarateC4(28*24:29*24))
% errorbar(Measured_hour,Measured_OA,sd_OA,'o')
% axis([0 24 0 50])
% title('Malate + Fumarate profile')
% xlabel('Hour')
% ylabel('Malate + Fumarate level (micromol C4 per g FW)')
% hold off
% 
% figure
% plot(0:24,sugar(28*24:29*24))
% title('Sucrose profile')
% xlabel('Hour')
% ylabel('Sucrose level (micromol per g FW)')
% %axis([0 24 0 6])
% 
% figure
% plot(1/24:1/24:40,sugar)
% title('Sucrose profile')
% xlabel('DAS')
% ylabel('Sucrose level (micromol per g FW)')
% 
% figure
% hold on
% errorbar(Measured_day,Measured_biomass,Measured_sd_biomass,'o')
% plot(1/24:1/24:40,FW)
% title('Rosette FW')
% xlabel('DAS')
% ylabel('g FW')
% hold off

end

FW'
starchC6(28*24:29*24)'
malatefumarateC4(28*24:29*24)'