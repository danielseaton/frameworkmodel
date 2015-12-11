function L = light_conditions(t,c)

t = mod(t, c.period);

LightOffset = 0; %Shifts light function up or down.
twilightPer = 0.05; %The duration of time between value of force in dark and value of force in light.
LightAmp = 1; %The amplitude of the light wave.

if (c.photoperiod == 0)
    L = 0;
elseif (c.photoperiod == c.period)
    L = 1;
else
    L = LightOffset + 0.5*LightAmp*(1 + tanh((c.period/twilightPer)*((t+c.dawn)/c.period - floor(floor(t+c.dawn)/c.period)))) - 0.5*LightAmp*(1 + tanh((c.period/twilightPer)*((t+c.dawn)/c.period - floor(floor(t+c.dawn)/c.period)) - c.photoperiod/twilightPer)) + 0.5*LightAmp*(1 + tanh((c.period/twilightPer)*((t+c.dawn)/c.period - floor(floor(t+c.dawn)/c.period)) - c.period/twilightPer));
end