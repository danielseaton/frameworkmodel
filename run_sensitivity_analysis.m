temp = 20.5;
rise = 0;
set = 12;
co2 = 42;
light = 145;


load('parameter.mat')

global p

p=parameter;
p0 = p;

deltaP = 0.01;

nP = length(p);
% nP = 2;

clock_genotype = {''};
% clock_genotype = {'prr9','prr7'};
clock_parameters = P2011_parameter_call(clock_genotype);


output_basal = integrationrebalMF_clock_for_sens_func(temp,rise,set,co2,light,clock_parameters);

nD = length(output_basal); %number of dimensions of output

sens = zeros(nP,nD);
errors = [];
for i = 1:nP
    i
    p = p0;
    p(i) = p(i)*(1+deltaP);
    try
        output = integrationrebalMF_clock_for_sens_func(temp,rise,set,co2,light,clock_parameters);
        sens(i,:) = ((output-output_basal)/deltaP)./output;
    catch
        errors = [errors;i];
    end
end

save('sens_analysis_results','output_basal','temp','rise','set','co2','light','p','deltaP','errors')

