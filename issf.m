function [ data ] = issf( time, minimum, maximum, on, off, twilight)
% Calculates 24 hour periodic data uses to simulate environmental inputs
%   Specification taken from http://jbr.sagepub.com/content/27/4/328.full

    t = mod(time-on,24);
    Tp = off-on;
    Tc = 24;
    T = max(0.0001,twilight);
    data = minimum + 0.5 * (maximum-minimum) * ( ( 1 + tanh(t/T)) - (1+ tanh((t-Tp)/T)) + (1 + tanh((t-Tc)/T)));

end
