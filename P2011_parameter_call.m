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