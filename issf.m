function [ data ] = issf( time, day, night, dayLength, twilight)
% Calculates 24 hour periodic data uses to simulate environmental inputs
%   Specification taken from http://jbr.sagepub.com/content/27/4/328.full

    % Move the center of the transition down by 0.5 hours to allow going from high to low in one hour
    % Also remove 1 because the simulation starts at time 1 but modular arithmetic has base 0 
    tPrime = mod(time-1+0.5,24);
    theta0 = night;
    theta1 = day - night;
    Tp = dayLength;
    Tc = 24;
    T = max(0.0001,twilight);
    
    data = theta0 + 0.5 * theta1 * ( ( 1 + tanh(tPrime/T)) - (1+ tanh((tPrime-Tp)/T)) + (1 + tanh((tPrime-Tc)/T)));

end
