function [S_intercept,Func_area] = calculate_area(t,S,petiole,input_size,Appear,Leaf_no,Sf,S_intercept,Func_area)
%% calculate area at a previous time point
%
% Input:
%   t - current time
%   S - leaf areas over time - n_timesteps x n_leaves
%   petiole - petiole factor determining crowding
%   input_size - number of timesteps of simulation
%   Appear - time of appearance of each leaf
%   Leaf_no - number of leaves over time
%   Sf - functioning leaf area over time - n_timesteps x n_leaves
%   Func_area - total functional leaf area over time
%
% Output:
%   S_intercept - total leaf area intercepting light
%   Func_area - total functional leaf area over time, updated within this function

if  Leaf_no(t-1) <= 15

    S_intercept(t-1) = sum(S(t-1,1:Leaf_no(t-1))); %Projected rosette area

else
    maxblade = sort(S(t-1,:),'descend');          
    crowding = (petiole-1)/(input_size-Appear(15))*(t-Appear(15))+1;%to decrease crowding at increasing LN
    S_intercept(t-1) = crowding*sum(maxblade(1:13));                              
end      


%Total functioning area (at the previous time point):

Func_area(t-1)=sum(Sf(t-1,1:Leaf_no(t-1)));
