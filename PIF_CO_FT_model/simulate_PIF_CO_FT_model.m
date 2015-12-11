function [dailyFTarea,FT_module_state] = simulate_PIF_CO_FT_model(sunrise,sunset,clock_output,clock_genotype,FT_module_state)

temperature = 22; %default parameter settings - may only take values of 22 or 27
flowering_genotype = {''}; %use to set e.g. CO or FKF1 mutants

% Load light conditions into 'c' for common light function
c.period = 24;
c.phase = 0;
c.dawn = sunrise;
c.photoperiod = sunset-sunrise;

% % Include model folders in path
% addpath('P2011_model')
% addpath('PIF_CO_FT_model')

parameters.clock = load_P2011_parameters(clock_genotype);

% Convert clock dynamics into input struct for PIF_CO_FT model
u = wrap_P2011_model_dynamics(clock_output.T,clock_output.Y,parameters.clock);

% Load PIF_CO_FT parameters
parameters.PIF_CO_FT = load_PIF_CO_FT_parameters(flowering_genotype,temperature);


% Run CO-PIF-FT model
% Simulate for one day
[T,Y] = ode15s(@(t,y) PIF_CO_FT_dynamics(t,y,parameters.PIF_CO_FT,u,c),[0 c.period],FT_module_state);
FT_module_state = Y(end,:)';


dailyFTarea = trapz(T,Y(:,15));