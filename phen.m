function [DayPhenThrm,FT_module_state] = phen(T,sunrise,sunset,geno,clock_parameters,clock_output,FT_module_state,p)

% Fixed parameters
%_________________

Tb = p(1);
Tvmin = p(2);
Tvmax = p(3);
Vsat = p(4);
v = p(5);
sigma = p(6);
m = p(7);
Dld = p(8);
CSDL = p(9);
CLDL = p(10);

if  geno ==1
    
    Fb = p(11);
    Dsd = p(12);
    Threshold = p(13);%2907;%3259
    Night = p(14);
    Phot_a = p(73);
    Phot_b = p(74);
    Phot_c = p(75);
    Phot_n = p(76);
    
else

    Fb = p(15);
    Dsd = p(16);
    Threshold = p(17);%3212;%4218
    Night = p(18);
    Phot_a = p(77);
    Phot_b = p(78);
    Phot_c = p(79);
    Phot_n = p(80);

end

%Overwrite photoperiod response parameters
Phot_a = 1;
Phot_b = -0.3937;
Phot_c = 2.8894;
Phot_n = 4;

%Calculation begins
%__________________
Thrm = zeros(24,1);

[dailyFTarea,FT_module_state] = simulate_PIF_CO_FT_model(sunrise,sunset,clock_output,clock_parameters,FT_module_state);

hour = 1:24;

for i = 1:24

    %Calculating thermal component
    %_____________________________

    if      sunrise>=hour(i) || sunset<=hour(i)-1
            fraction_light(i)=0;
    elseif  sunrise<=hour(i)-1 && sunset>hour(i)
            fraction_light(i)=1;
    elseif  sunrise>=hour(i)-1
            fraction_light(i)=hour(i)-sunrise;
    else
            fraction_light(i)=sunset-hour(i)+1;
    end



    if      fraction_light(i)==0
            Thrm(i)=Night*max(T(i)-Tb,0);
    elseif  fraction_light(i)==1    
            Thrm(i)=max(T(i)-Tb,0);
    else 
            Thrm(i)=max(0,(T(i)-Tb)*fraction_light(i)) + Night*max(0,(T(i)-Tb)*(1-fraction_light(i)));
    end




    %Calculating photoperiod component
    %_________________________________

%     dl(i) = sunset(i) - sunrise(i);

%     if  i == 1
% 
%         [FTarea24,yo] = link(dl(1),sunrise(1)); %initialise
%         [FTarea1,yo] = sublink(hour(i),dl(i),sunrise(i),yo);
%         FTarea24 = [FTarea24(2:end) FTarea1];
%         dailyFTarea(i) = sum(FTarea24);                                       
%     else
%         [FTarea1,yo] = sublink(hour(i),dl(i),sunrise(i),yo);
%         FTarea24 = [FTarea24(2:end) FTarea1];
%         dailyFTarea(i) = sum(FTarea24);
%     end

%     dailyFTarea(i) = 0;
    Phot(i) = Phot_a + Phot_b*Phot_c^Phot_n/(Phot_c^Phot_n+dailyFTarea^Phot_n);
%     Phot(i)



    %Calculating vernalization effectiveness
    %_______________________________________


    if	T(i)>=Tvmin && T(i)<=Tvmax

        Ve(i)=exp(v)*(T(i)-Tvmin)^m*(Tvmax-T(i))^sigma*1;
    else	
        Ve(i)=0;
    end



    %Calculating cumulative vernalization hours
    %__________________________________________


    Vh(i)=sum(Ve(1:i));




    %Calculating vernalization fraction
    %__________________________________

    if	Vh(i)<=Vsat

        Vern(i) = Fb + Vh(i)*(1 - Fb)/Vsat;
    else
        Vern(i)= 1;	
    end



    %Calculating modified photothermal unit
    %______________________________________

    mptu(i) = Vern(i)*Phot(i)*Thrm(i);

end


%Output:
%-------
DayPhenThrm=sum(mptu);