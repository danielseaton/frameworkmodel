function [] = run_simulation(resultFileName,tempDay,tempNight,tempDayLength,tempTwilight,co2Day,co2Night,co2DayLength,co2Twilight,lightDay,lightNight,lightDayLength,lightTwilight)

  fileID = fopen(resultFileName,'w');

  %hour,T,sunrise,sunset,CO2,PAR
  load('weather.mat')

  hour=weather(:,1);

  rise = 1;
  set = rise+lightDayLength;

  T=issf(hour, tempDay, tempNight, tempDayLength, tempTwilight);
  sunrise=rise*weather(:,3);
  sunset=set*weather(:,4); 
  CO2=issf(hour, co2Day, co2Night, co2DayLength, co2Twilight);
  PAR=issf(hour, lightDay, lightNight, lightDayLength, lightTwilight);
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

  d = 1-w; %dry matter

  load('parameter.mat')

  clock_genotype = {''};
  % Could specify mutants here, e.g.: clock_genotype = {'prr9','prr7'};
  clock_parameters = P2011_parameter_call(clock_genotype);

  starch_genotype = {''};
  starch_parameters = starch_parameter_call(starch_genotype);


  simulate_FM(hour,T,sunrise,sunset,CO2,PAR,Photoperiod,clock_parameters,starch_parameters,parameter,d,mf_use,true,fileID);

  % Close the result file
  fclose(fileID);

end
