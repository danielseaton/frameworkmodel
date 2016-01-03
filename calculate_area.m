function [S_intercept,Func_area,Destructive_projected_area] = calculate_area(t,S,petiole,input_size,Appear,Leaf_no,Sf)

%FUNC: calculate area a previous time point
if  Leaf_no(t-1) <= 15

    S_intercept(t-1) = sum(S(t-1,1:Leaf_no(t-1))); %Projected rosette area

else
    maxblade = sort(S(t-1,:),'descend');          
    crowding = (petiole-1)/(input_size-Appear(15))*(t-Appear(15))+1;%to decrease crowding at increasing LN
    S_intercept(t-1) = crowding*sum(maxblade(1:13));                              
end      


%Total functioning area (at the previous time point):

Func_area(t-1)=sum(Sf(t-1,1:Leaf_no(t-1))); 
Destructive_projected_area(t-1) = sum(S(t-1,1:Leaf_no(t-1)));