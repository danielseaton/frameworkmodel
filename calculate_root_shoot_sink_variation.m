function [fL,fR,rsratio] = calculate_root_shoot_sink_variation(t,Leaf_no,TLcot,Maxfcot,...
    TLtrue,MaxfL,aL,bL,CumThrm,Appear,TR,MaxfR,aR,bR,PR,PL,rsratio,leaf_factor,root_factor,fL,fR)

%To calculate sink variation for leaf:
for i = 1:Leaf_no(t)

    if i <= 2
       TL = TLcot;
       Maxf = Maxfcot;

       fL(t,i) = ((CumThrm(t)-CumThrm(Appear(i))+0.5)/TL)^(aL-1)...
              *(1-(CumThrm(t)-CumThrm(Appear(i))+0.5)/TL)^(bL-1);

       if real(fL(t,i))~=fL(t,i);
          fL(t,i) = 0; 
       end      

       fL(t,i) = fL(t,i)/Maxf;        

    else
       TL =TLtrue;            
       Maxf = MaxfL;

       fL(t,i) = ((CumThrm(t)-CumThrm(Appear(i))+0.5)/TL)^(aL-1)...
              *(1-(CumThrm(t)-CumThrm(Appear(i))+0.5)/TL)^(bL-1);

       if real(fL(t,i))~=fL(t,i);
          fL(t,i) = 0; 
       end      

       fL(t,i) = fL(t,i)/Maxf;

    end


end


%To calculate sink variation for root:       
fR(t) = ((CumThrm(t)+0.5)/TR)^(aR-1)*(1-(CumThrm(t)+0.5)/TR)^(bR-1);

fR(t) = fR(t)/MaxfR; %sink strength   


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

num   = fR(t)*PR*rootconversion;
denom = sum(fL(t,:))*PL*leafconversion;
rsratio(t) = num/denom;