function [net_rate] = photosynthesis(CO2,Tleaf,PAR,vlmax25,daylength,p)
%Calculate rate of photosynthesis limited by Rubisco
%
% Inputs:
%    CO2 - CO2 partial pressure (Pa)
%    Tleaf - Leaf temperature (degrees Centigrade)
%    PAR - Total absorbed photosynthetically active radiation per unit leaf area (micromol m-2 s-1)
%    vlmax25 - Initial photosynthetic rubisco capacity per unit leaf area at 25degC (micromol CO2 m-2 s-1)
%    daylength = sunset - sunrise, in h
%    p - Vector of parameters
%
% Outputs:
%    net_rate - Net rate of leaf photosynthesis (micromol CO2 m-2 leaf s-1)
%
% Functions called:
%    none
%
%   Copyright 2018 Yin Hoon Chew, Daniel Seaton, Andrew Millar, and The University of Edinburgh
%
%   Licensed under the Apache License, Version 2.0 (the "License");
%   you may not use this file except in compliance with the License.
%   You may obtain a copy of the License at
%
%       http://www.apache.org/licenses/LICENSE-2.0
%
%   Unless required by applicable law or agreed to in writing, software
%   distributed under the License is distributed on an "AS IS" BASIS,
%   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%   See the License for the specific language governing permissions and
%   limitations under the License.


%Rasse and Tocquin, 2006


%___________________________________________________


%Rubisco activity for CO2
%________________________

R = p(39); %gas constant
act_ener_kc25 = p(40); %activation energy for kc25 (kJ mol-1)
KmRubisco25 = p(41); %Michaelis-Menten constant of Rubisco for CO2 at 25(Pa)

numC = act_ener_kc25*(Tleaf-25);
denomC = 298*R*(Tleaf+273);

kc_Tleaf = KmRubisco25*exp(numC/denomC); %M-M Rub co2 at leaf T


%Rubisco activity for O2
%_______________________

act_ener_ko25 = p(42); %activation energy for ko25
KmRubO = p(43); %Michaelis-Menten constant of Rubisco for O2 at 25 (Pa) 

numO = act_ener_ko25*(Tleaf-25);

ko_Tleaf = KmRubO*exp(numO/denomC); %M-M constant for O2 at leaf T


%Photosynthetic Rubisco capacity at leaf T
%_________________________________________

act_ener_vlmax = p(44); %activation energy for vlmax25

numR = act_ener_vlmax*(Tleaf-25);

vlmax_Tleaf = vlmax25*exp(numR/denomC); %micromol CO2 m-2 s-1


%Effective M-M constant of Rubisco (Pa)
%______________________________________

o2_parpre = p(45); %O2 partial pressure (Pa)

kprim_effective = kc_Tleaf*(1+o2_parpre/ko_Tleaf);


%Rubisco-limited rate
%____________________

co2_comp_point = p(46)+p(47)*(Tleaf-25)+p(48)*(Tleaf-25)^2; %CO2 compensation point in the absence of mitochondrial respiration (Pa)
int_co2 = p(49)*CO2; %intercell CO2 partial pressure (Pa)

numAV = int_co2-co2_comp_point;
denomAV = int_co2+kprim_effective;

av_rub = vlmax_Tleaf*(numAV/denomAV);



%Calculate rate of photosynthesis limited by RuBP regeneration
%_____________________________________________________________
%_____________________________________________________________


%Electron transport
%__________________

act_ener_jm25 = p(50); %activation energy for jm25 (kJ mol-1)
h_curvature = p(51); %Curvature parameter of Jm (J mol-1)
s_elec = p(52); %Electron transport temperature response parameter (J K-1 mol-1)

numt1 = act_ener_jm25*(Tleaf-25);
denomt1 = 298*R*(Tleaf+273);
term1 = exp(numt1/denomt1);

numt2 = s_elec*298-h_curvature;
denomt2 = R*298;
term2 = 1+exp(numt2/denomt2); %intermediate calculation for jm

numt3 = s_elec*(273+Tleaf)-h_curvature;
denomt3 = R*(273+Tleaf);
term3 = 1+exp(numt3/denomt3); %intermdeiate calculation fro jm

daylengths = [4,6,8,12,14,18];
jv_values = [2.1,2.1,2.1,1.7,1.4,1.2];
%linear interpolation to calculate jv as a function of daylength (in the
%range 4-18hrs)
jv = interp1(daylengths,jv_values,daylength);     

jm25_pot = jv*vlmax25;

jm_pot = jm25_pot*term1*term2/term3; %Potential rate of electron transport oer unit leaf area (micromol m-2 s-1)


%Photosystem II
%______________

fspec = p(54); %Spectral correction factor. PAR absorbed by tissue other than the chloroplast lamella

ile_par = PAR*(1-fspec)/2; %PAR effectively absorbed by PSII

cc = ile_par*jm_pot;
bb = -1*(ile_par+jm_pot);


%jl electron transport (quadratic solution)
%__________________________________________

thetal_curve = p(55); %Curvature of leaf response of electron transport to irradiance
aa = thetal_curve;

rho = bb^2-4*aa*cc;

if  rho <= 0
    sol = 0;
else
    sol = (-1*bb-rho^0.5)/(2*aa); %Solving the quadratic equation thetal*j^2 - (ile+jm)*j + ile*jm = 0
end


%Rubisco regeneration
%____________________

numAJ = sol*(int_co2-co2_comp_point);
denomAJ = 4*(int_co2+2*co2_comp_point);

aj_RuBP = numAJ/denomAJ; %Rate of phtoosynthesis limited by RuBP
    

if  rho <= 0
    net_rate = av_rub;
else
    net_rate = min(av_rub,aj_RuBP); %Net rate of leaf photosynthesis (micromol CO2 m-2 leaf s-1)
end

