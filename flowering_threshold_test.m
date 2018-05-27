function has_flowered = flowering_threshold_test(CumPhenThrm,geno,p)
%% determine whether the flowering threshold has been passed
%
% Input:
%   CumPhenThrm - cumulative thermal units, adjusted by the phenology model depending on the photoperiod
%   geno - genotype specified. Different values result in different thresholds for flowering (for Col or Ler)
%   p - vector of FM parameters
%
% Output:
%   has_flowered - boolean, indicating whether the flowering threshold has been passed



if  geno ==1
    
    Threshold = p(13);%2907;%3259
    
else
    Threshold = p(17);%3212;%4218

end

if CumPhenThrm > Threshold
    has_flowered = true;
else
    has_flowered = false;
end
