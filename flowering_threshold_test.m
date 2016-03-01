function has_flowered = flowering_threshold_test(CumPhenThrm,geno,p)

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