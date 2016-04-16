temp = 20.5;
rise = 0;
set = 12;
co2 = 42;
light = 145;

%Calling for meteorological data
%_______________________________  

%hour,T,sunrise,sunset,CO2,PAR
load('weather.mat')

hour=weather(:,1);
T=temp*weather(:,2);
sunrise=rise*weather(:,3);
sunset=set*weather(:,4); 
CO2=co2*weather(:,5); %CO2 partial pressure (Pa)
PAR=light*weather(:,6); %total absorbed PAR per unit leaf area (micromol m-2 s-1)
Photoperiod = sunset(end) - sunrise(end);



%Specifying the genotype for the clock and starch models
experiment_genotype_index = 1;
switch experiment_genotype_index
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


load('parameter.mat')

p=parameter;

clock_genotype = {''};
% clock_genotype = {'prr9','prr7'};
starch_genotype = {''};

clock_parameters = P2011_parameter_call(clock_genotype);
nCP = length(clock_parameters);

starch_parameters = starch_parameter_call(starch_genotype);
nSP = length(starch_parameters);

[output,sim_data] = run_simulation(hour,T,sunrise,sunset,CO2,PAR,Photoperiod,clock_parameters,starch_parameters,p0,d,mf_use);
