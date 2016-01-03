function [Appear,Leaf_mass,Si,Leaf_no,GC,GCstart] = determine_current_leaf_number(t,probability,CumThrm,phyllochron,Appear,Leaf_mass,Si,Leaf_no,GC,GCstart)

if  (CumThrm(t) - CumThrm(GCstart(end))) >= phyllochron % A growth cycle has been reached

    GC = GC+1;
    GCstart(GC+1) = t; % the timepoint when a new growth cycle starts

    Leaf_no(t) = Leaf_no(t-1) + binornd(1,probability); % A new leaf may/may not be initiated

    if Leaf_no(t) ~= Leaf_no(t-1) % A new leaf is initiated
       Appear(Leaf_no(t)) = t;
       Leaf_mass(t-1,Leaf_no(t))=0;
       Si(t-1,Leaf_no(t))=0;
    end              

else
    Leaf_no(t) = Leaf_no(t-1);        
end          
