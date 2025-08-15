% Updated 01/12/2025
% For citations regarding vector parameters, please refer to comments in West Nile code written by Dr. Wandi Ding and Dr. Rachel Leander.

function [pU] = model_DE_Parameters_Control(larvicide_type, Tf)

% Parameter values for simulations:
pU = zeros(1,12);
pV=model_DE_Parameters_Vectors();


% INSECTICIDE PARAMETERS

% 1 CORRESPONDS TO METHOPRENE
if larvicide_type == 1
%%%%%%%Computation of p(46) and p(47) for s-methoprene briquets%%%%%%%%%%%%

min_ef=.03; % Rough approximation of min_ef and min_ef days

% The product assessment has the product lasting up to 150 days and 69.5% effective at 120 days
% It does not last this long in the field. 
min_ef_day=150;
half_ef_day=100;
max_ef=1/(1+(min_ef/(1-min_ef))^(half_ef_day/(min_ef_day-half_ef_day)));

pU(1) = max_ef*(pV(8)+pV(7))/(1-max_ef); % max rate at which larvicide kills larvae
pU(2)= -log((1-max_ef)/max_ef)/half_ef_day; % rate at which larvicide decays

end

% 2 CORRESPONDS TO VECTOBAC
if larvicide_type == 2
%%%%%%%Computation of p(46) and p(47) for vectobac%%%%%%%

min_ef = .22;
min_ef_day = 42;
mid_ef_day = 34;
mid_ef = .57;

max_ef=1/(1+((min_ef/(1-min_ef))^(mid_ef_day/(min_ef_day-mid_ef_day)))*((1-mid_ef)/mid_ef)^(min_ef_day/(min_ef_day-mid_ef_day)));
pU(1) = max_ef*(pV(8)+pV(7))/(1-max_ef); % max rate at which larvicide kills larvae
pU(2)= -log(mid_ef*(1-max_ef)/((1-mid_ef)*max_ef))/mid_ef_day; % rate at which larvicide decays

end

pU(3) = 24; % adulticide decay rate
pU(4) = -log(0.1)*pU(3)*2; % max rate at which adulticide kills adult vectors
per_remain_one_hour = exp(pU(3)*0.5*exp(-pU(4)/24)/pU(4)-pU(3)*0.5/pU(4)); % Percent of adulticide remaining after one hour

% CONTROL PARAMETERS

pU(5) = 5000;   % the weight of the cost of the infected components or vectors,
                % according to disease or vector respectively, in the integral cost term.

pU(6) = 1; % weight of cost of larvacide
pU(7) = 10; % weight of cost of adulticide

pU(8) = 0.05; % weight of cost of time
pU(9) = 5000; % cost of eggs at the final time
pU(10) = -100000; % cost of hosts at the final time; value used in paper simulations

pU(11) = Tf; % maximum time between controls
pU(12) = 1; % minimum time between controls

end