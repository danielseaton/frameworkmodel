function [clock_output] = circadian_module(sunrise,sunset,clock_state_0,clock_genotype)

%clock_parameters

% Clock parameters:
P = P2011_parameter_call(clock_genotype);

Tout = [0:0.05:27];
% Run model for 27 hours:
[T,Y]=ode15s(@(t,y) P2011_dynamics(t,y,P,sunrise,sunset),Tout,clock_state_0);

clock_output = struct();
clock_output.T = T;
clock_output.Y = Y;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P2011 clock model dynamics
function Func = P2011_dynamics(t,y,P,sunrise,sunset)

    P = num2cell(P);
    [q1,q2,q3,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,...
        p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29,p30,p31,m1,m2,m3,m4,m5,m6,...
        m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24,...
        m25,m26,m27,m28,m29,m30,m31,m32,m33,m34,m35,m36,m37,m38,m39,n1,n2,n3,...
        n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,...
        g11,g12,g13,g14,g15,g16,a,b,c,d,e,f] = deal(P{:});

    Func = zeros(30, 1);

    cond.period=24;
    cond.dawn=sunrise;
    cond.photoperiod = sunset-sunrise;
    L = light_conditions(t,cond);
%     LightOffset = 0; %Shifts light function up or down.
%     twilightPer = 0.05; %The duration of time between value of force in dark and value of force in light.
%     LightAmp = 1; %The amplitude of the light wave.
%     period = 24;
% 
%     t = mod(t,period);
% 
%     L = LightOffset + 0.5*LightAmp*(1 + tanh((period/twilightPer)*((t+sunrise)/period - ...
%         floor(floor(t+sunrise)/period)))) - 0.5*LightAmp*(1 + ...
%         tanh((period/twilightPer)*((t+sunrise)/period - floor(floor(t+sunrise)/period)) - ...
%         (sunset-sunrise)/twilightPer)) + 0.5*LightAmp*(1 + ...
%         tanh((period/twilightPer)*((t+sunrise)/period - floor(floor(t+sunrise)/period)) -...
%         period/twilightPer));
    
    %ODEs
    % y(1)  LHY mRNA                   
    % y(2)  P
    % y(3)  GI-ZTL 
    % y(4)  GI-ELF3 cytoplasm
    % y(5)  LHY prot
    % y(6)  TOC1 mRNA
    % y(7)  PRR9 prot
    % y(8)  PRR5 (NI) mRNA
    % y(9)  PRR5 (NI) prot
    % y(10) GI prot cytoplasm
    % y(11) TOC1 prot
    % y(12) ZTL
    % y(13) EC
    % y(14) GI mRNA
    % y(15) PRR9 mRNA
    % y(16) PRR7 mRNA
    % y(17) PRR7 prot
    % y(18) ELF4 mRNA
    % y(19) ELF4 prot
    % y(20) LHY prot modif.
    % y(21)     HY5 prot used only for fitting of COP1 parameters
    % y(22)     HFR1 prot used only for fitting of COP1 parameters
    % y(23) ELF3 mRNA
    % y(24) ELF3 cytoplasm
    % y(25) ELF3 nuclear
    % y(26) COP1 nuclear night
    % y(27) COP1 nuclear day
    % y(28) LUX mRNA
    % y(29) LUX prot
    % y(30) COP1 cytoplasm

    Gn=p28*y(10)/(p29+m19+p17*y(25));
    EGn=(p18*y(4)+p17*y(25)*Gn)/(m9*y(26)+m10*y(27)+p31);
    e34=p25*y(19)*y(25)/(p26*y(29)+p21+m36*y(26)+m37*y(27));

    Func(1) = g1^a/(g1^a+(y(7)+y(17)+y(9)+y(11))^a)*(q1*L*y(2)+n1)-y(1)*(m1*L+m2*(1-L));
    Func(2) = p7*(1-L)*(1-y(2))-m11*y(2)*L;
    Func(3) = p12*L*y(12)*y(10)-p13*y(3)*(1-L)-m21*y(3);
    Func(4) = p17*y(24)*y(10)-m9*y(4)*y(30)-p18*y(4)+p31*EGn;
    Func(5) = (p2+p1*L)*y(1)-m3*y(5)-p3*y(5)^c/(y(5)^c+g3^c);
    Func(6) = 1*n2*g4/(g4+y(13))*(g5^e/(g5^e+y(5)^e))-y(6)*m5;
    Func(7) = p8*y(15)-(m13+m22*(1-L))*y(7);
    Func(8) = 1*(n10*y(20)^e/(g12^e+y(20)^e)+n11*y(17)^b/(g13^b+y(17)^b))-m16*y(8);
    Func(9) = p10*y(8)-(m17+m24*(1-L))*y(9);
    Func(10)= p11*y(14)-m19*y(10)-p12*L*y(12)*y(10)+p13*y(3)*(1-L)-p17*y(24)*y(10)-p28*y(10)+p29*Gn;
    Func(11)= p4*y(6)-m8*y(11)-(m6+m7*(1-L))*y(11)*(p5*y(12)+y(3));
    Func(12)= 1*p14-m20*y(12)-p12*L*y(12)*y(10)+p13*y(3)*(1-L);
    Func(13)= p26*y(29)*e34-m36*y(13)*y(26)-m37*y(13)*y(27)-m32*y(13)*(1+p24*L*(EGn+Gn)^d/(g7^d+(EGn+Gn)^d));
    Func(14)= 1*(q2*L*y(2)+g15^e/(g15^e+y(5)^e)*g14/(g14+y(13))*n12)-y(14)*m18;
    Func(15)= 1*q3*L*y(2)+g8/(g8+y(13))*(n4+n7*y(5)^e/(y(5)^e+g9^e))-m12*y(15);
    Func(16)= 1*(n8*(y(5)+y(20))^e/(g10^e+(y(5)+y(20))^e)+n9*y(7)^f/(g11^f+y(7)^f))-m14*y(16);
    Func(17)= p9*y(16)-y(17)*(m15+m23*(1-L));
    Func(18)= n13*(g6^e/(g6^e+y(5)^e))*g2/(g2+y(13))-y(18)*m34;
    Func(19)= p23*y(18)-m35*y(19)-p25*y(25)*y(19)+p21*e34;
    Func(20)= p3*y(5)^c/(y(5)^c+g3^c)-m4*y(20);
    Func(21)= p22-m38*y(21)*y(27)-m25*y(21)*y(26);
    Func(22)= p30-m28*y(22)*y(26);
    Func(23)= 1*n3*g16^e/(g16^e+y(5)^e)-m26*y(23);
    Func(24)= p16*y(23)-m9*y(24)*y(30)-p17*y(24)*y(10)-p19*y(24)+p20*y(25);
    Func(25)= p19*y(24)-p20*y(25)-m29*y(25)*y(26)-m30*y(25)*y(27)-p25*y(25)*y(19)+p21*e34-p17*y(25)*Gn;
    Func(26)= p6*y(30)-n6*L*y(2)*y(26)-n14*y(26)-m27*y(26)*(1+p15*L);
    Func(27)= 1*(n14*y(26)+n6*L*y(2)*y(26))-m31*(1+m33*(1-L))*y(27);
    Func(28)= n13*(g6^e/(g6^e+y(5)^e))*g2/(g2+y(13))-y(28)*m34;
    Func(29)= p27*y(28)-m39*y(29)-p26*y(29)*e34;
    Func(30)= 1*n5-p6*y(30)-m27*y(30)*(1+p15*L);
end

function P = P2011_parameter_call(mutant)
% this function returns parameters for P2011 according to the desired
% mutant.

% "mutant" is a cell with entries for each gene desired to be knocked out

% basal parameter set:
m1=0.54; m2=0.24; m3=0.2; m4=0.2; m5=0.3; m6=0.3; m7=0.7; m8=0.4; m9=1.1; m10=1; m11=1; m12=1; m13=0.32; m14=0.4; m15=0.7; 
m16=0.5; m17=0.5; m18=3.4; m19=0.2; m20=0.6; m21=0.08; m22=0.1; m23=1.8; m24=0.1; m25=1.8; m26=0.5; m27=0.1; m28=20; m29=5; 
m30=3; m31=0.3; m32=0.2; m33=13; m34=0.6; m35=0.3; m36=0.1; m37=0.8; m38=0.5; m39=0.3; 

n1=2.6; n2=0.64; n3=0.29; n4=0.07; n5=0.23; n6=20; n7=0.2; n8=0.5; n9=0.2; n10=0.4; n11=0.6; n12=12.5; n13=1.3; n14=0.1; n15=0.4; n16=0.23; n17=0.1; 

p1=0.13; p2=0.27; p3=0.1; p4=0.56; p5=4; p6=0.6; p7=0.3; p8=0.6; p9=0.8; p10=0.54; p11=0.51; p12=3.4; p13=0.1; p14=0.14; p15=3;
p16=0.62; p17=4.8; p18=4; p19=1; p20=0.1; p21=1; p22=0.5; p23=0.37; p24=10; p25=8; p26=0.3; p27=0.8; p28=2; p29=0.1; p30=0.9; p31=0.1;

g1=0.1; g2=0.01; g3=0.6; g4=0.01; g5=0.15; g6=0.3; g7=0.6; g8=0.01; g9=0.3; g10=0.5; g11=0.7; g12=0.2; g13=1; g14=0.004; g15=0.4;
g16=0.3; g17=0.03; g18=0.3; g19=0.5; g20=1.1; g21=0.5; g22=1.5; g23=0.01; g24=3; g25=0.03; g26=0.3; g27=0.1; g28=0;

a=2; b=2; c=2; d=2; e=2; f=2; g=0; h=0; 

q1=1.2; q2=1.56; q3=2.8; q4=0; 

% modify according to mutant:

for i = 1:length(mutant)
    if strcmp(mutant{i},'lhy')
        % LHY/CCA1 mutant{i}:
        q1 = 0;
        n1 = 0;
    elseif strcmp(mutant{i},'CCA1OX')
        % CCA1-OX mutant{i}:
        q1 = 0;
        n1 = 5;
        g1 = 10000;
    elseif strcmp(mutant{i},'gi')
        % GI mutant{i}:
        q2 = 0;
        n12 = 0;
    elseif strcmp(mutant{i},'prr5')
    %   PRR5 mutant{i}:
        n10 = 0;
        n11 = 0;
    elseif strcmp(mutant{i},'prr7')
        %PRR7 mutant{i}:
        n8 = 0;
        n9 = 0;
    elseif strcmp(mutant{i},'prr9')
        %PRR9 mutant{i}:
        q3 = 0;
        n4 = 0;
        n7 = 0;
    elseif strcmp(mutant{i},'toc1')
        %TOC1 mutant{i}:
        n2 = 0;
    elseif strcmp(mutant{i},'elf3')
        %ELF3 mutant{i}:
        n3 = 0; 
    elseif strcmp(mutant{i},'ztl')
        %ZTL mutant{i}:
        p14 = 0;
    end
end

P = [q1,q2,q3,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18...
    p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29,p30,p31,m1,m2,m3,m4,m5,m6...
    m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23,m24...
    m25,m26,m27,m28,m29,m30,m31,m32,m33,m34,m35,m36,m37,m38,m39,n1,n2,n3...
    n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,g1,g2,g3,g4,g5,g6,g7,g8,g9,g10...
    g11,g12,g13,g14,g15,g16,a,b,c,d,e,f];
end


end