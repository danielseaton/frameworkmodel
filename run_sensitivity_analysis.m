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

output_filename = input(['Enter name for datafile output \n(e.g. sens_analysis_results_WT):'], 's');


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
p0 = p;

deltaP = 0.01;

nP = length(p);
% nP = 2;

clock_genotype = {''};
% clock_genotype = {'prr9','prr7'};
starch_genotype = {''};
clock_parameters0 = P2011_parameter_call(clock_genotype);
clock_parameters = clock_parameters0;
nCP = length(clock_parameters);

starch_parameters0 = starch_parameter_call(starch_genotype);
starch_parameters = starch_parameters0;
starch_parameter_names = fieldnames(starch_parameters);
nSP = length(starch_parameter_names);

run_phenology_model = 0;

[output_basal,~] = simulate_FM(hour,T,sunrise,sunset,CO2,PAR,Photoperiod,clock_parameters,starch_parameters,p0,d,mf_use,run_phenology_model);

nD = length(output_basal); %number of dimensions of output

sens = zeros(nP+nCP+nSP,nD);
errors = [];
for i = 1:nP
    i
    p = p0;
    p(i) = p(i)*(1+deltaP);
    try
        [output,~] = simulate_FM(hour,T,sunrise,sunset,CO2,PAR,Photoperiod,clock_parameters,starch_parameters,p,d,mf_use,run_phenology_model);
        sens(i,:) = ((output-output_basal)/deltaP)./output;
    catch
        errors = [errors;i];
    end
end

p = p0;
for i = nP+1:nP+nCP
    i
    clock_parameters = clock_parameters0;
    clock_parameters(i-nP) = clock_parameters(i-nP)*(1+deltaP);
    try
        [output,~] = simulate_FM(hour,T,sunrise,sunset,CO2,PAR,Photoperiod,clock_parameters,starch_parameters,p,d,mf_use,run_phenology_model);
        sens(i,:) = ((output-output_basal)/deltaP)./output;
    catch
        errors = [errors;i];
    end
end

p = p0;
clock_parameters = clock_parameters0;
for i = nP+nCP+1:nP+nCP+nSP
    i
    starch_parameters = starch_parameters0;
    starch_parameters.(starch_parameter_names{i-nP-nCP}) = starch_parameters.(starch_parameter_names{i-nP-nCP})*(1+deltaP);
    try
        [output,~] = simulate_FM(hour,T,sunrise,sunset,CO2,PAR,Photoperiod,clock_parameters,starch_parameters,p,d,mf_use,run_phenology_model);
        sens(i,:) = ((output-output_basal)/deltaP)./output;
    catch
        errors = [errors;i];
    end
end

save(output_filename,'output_basal','temp','rise','set','co2','light','p0','clock_parameters0','clock_genotype','deltaP','errors','sens')