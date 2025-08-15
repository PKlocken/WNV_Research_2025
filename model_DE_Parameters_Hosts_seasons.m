function hostParams = model_DE_Parameters_Hosts_seasons(seasons, urban, transmit, varargin)

% Inputs:
%   urban           - Urban scenario (1 = less urban, 2 = more urban, 3 = amplify, 4 = dilute, 5=make4 most preferred)
%   transmit        - Transmission toggles (e.g., 4, 34, 1, 0)
%   varargin        - Required if urban > 2 (deppends on urban value e.g., if urban =3,4 then 1 = change density, 2 = change preference)


if nargin < 2
    error('Urban and transmit inputs are required');
end


if urban > 2 && urban < 5 % Tag input for scenarios urban = 3 or 4
    if isempty(varargin)
        error('Additional input required for urban > 2 cases.');
    end
    tag_input = varargin{1};
    factor = 100;
elseif urban == 5
    if isempty(varargin)
        error('Additional input required for urban > 2 cases.');
    end
    factor = varargin{1};    
else
    tag_input = 0;
    factor = 1;
end

% Updated 02/21/2025
% For citations regarding vector parameters, please refer to comments in West Nile code written by Dr. Wandi Ding and Dr. Rachel Leander.

% Host Parameter values for simulations:

% Host type 1: vulnerable amplifiers with horizontal transmission
pH1 = zeros(1,10);
pH1(1) = 0.68; % host-to-mosquito transmission probability host group 1, p_hm1

pH1(3) = (1/3)*((1/4.7)+(log(0.6)/365))+log(0.6)/365; % WNV recovery rate host group 1, g1, gamma/(g+gamma+mu)=0.75
pH1(4) = (1/4.7)+(log(0.6)/365) ; % WNV induced death rate host group 1, gamma1, mu+gamma=1/4.7
%pH1(5) = (0.15-log(0.6))/365; % annual per capita birth rate host group 1, Lambda1
if seasons==2
    pH1(5) = 4*(0.15-log(0.6))/365; % 4 x annual per capita birth rate host group 1, Lambda1
    %assume all births occur during a three month breeding season
else
    pH1(5)=0;
end
pH1(6) = -log(0.6)/365; % natural death rate host group 1, mu_h1

pH1(2) = 3.5*(pH1(6)+pH1(4)+pH1(3));% maximal horizontal transmission rate host group 1, omega1,
%3.5=omega/(gamma+g+mu)
%pH1(7) = 0.0001; % Less Urban carrying capacity host group 1, c_h1 in m^2
%pH1(8) = 0.000025; % More Urban carrying capacity host group 1, c_h1 in
%m^2
pH1(7) = 1; % Less Urban carrying capacity host group 1, c_h1 in ha
pH1(8) = 0.25; % More Urban carrying capacity host group 1, c_h1 in ha
pH1(9) = 1; % biting preference host group 1, alpha1
pH1(10) = 1; % prob. horizontal transmission host group 1, p_hh1    %p_hh values all set to 1

% Host type 2: vulnerable amplifiers without horizontal transmission
pH2 = zeros(1,10);
pH2(1) = 0.53; % host-to-mosquito transmission probability host group 2, p_hm2
pH2(2) = 0; % maximal horizontal transmission rate host group 2, omega2
pH2(3) = (1/4.7)+2*(log(0.6)/365); % WNV recovery rate host group 2, g2, gamma/(g+gamma+mu)=0.5
pH2(4) = (1/4.7)+(log(0.6)/365); % WNV induced death rate host group 2, gamma2
%pH2(5) = (0.3-log(0.6))/365; % annual per capita birth rate host group 2, Lambda2
if seasons==2
    pH2(5) = 4*(0.3-log(0.6))/365;% 4 x annual per capita birth rate host group 1, Lambda1
%assume all births occur during a three month breeding season
else
    pH2(5)=0;
end
pH2(6) = -log(0.6)/365; % natural death rate host group 2, mu_h2
%pH2(7) = 0.0001; % Less Urban carrying capacity host group 2, c_h2
%pH2(8) = 0.0004; % More Urban carrying capacity host group 2, c_h2
pH2(7) = 2; % Less Urban carrying capacity host group 2, c_h2
pH2(8) = 5; % More Urban carrying capacity host group 2, c_h2
pH2(9) = 0.1; % biting preference host group 2, alpha2
%pH2(9) = 1; 
pH2(10) = 1; % prob. horizontal transmission host group 2, p_hh2 %p_hh values all set to 1

% Host type 3: invulnerable amplifier
pH3 = zeros(1,10);
pH3(1) = 0.36; % host-to-mosquito transmission probability host group 3, p_hm3
pH3(2) = 0; % maximal horizontal transmission rate host group 3, omega3
pH3(3) = 0.3333; % WNV recovery rate host group 3, g3
pH3(4) = 0; % WNV induced death rate host group 3, gamma3
%pH3(5) = (0.1570-log(0.6))/365; % per capita birth rate host group 3, Lambda3
if seasons==2
    pH3(5) = 4*(0.1570-log(0.6))/365; %4x annual birth rate, assume all births occur in a three month breeding season;
else
    pH3(5)=0;
end
pH3(6) = -log(0.6)/365; % natural death rate host group 3, mu_h3
%pH3(7) = 0.0002; % Less Urban carrying capacity host group 3, c_h3 (hosts
%per m^2)
%pH3(8) = 0.00005; % More Urban carrying capacity host group 3, c_h3 (hosts
%per m^2)
pH3(7) = 2; % Less Urban carrying capacity host group 3, c_h3 (hosts per ha)
pH3(8) = 1; % More Urban carrying capacity host group 3, c_h3 (hosts per ha)
pH3(9) = 1; % biting preference host group 3, alpha3
pH3(10) = 1; % host-host contact rate host group 3, p_hh3    %p_hh values all set to 1

% Host type 4: diluters
pH4 = zeros(1,10);
pH4(1) = 0.1; % host-to-mosquito transmission probability host group 4, p_hm4
pH4(2) = 0; % maximal horizontal transmission rate host group 4, omega4
pH4(3) = 1; % WNV recovery rate host group 4, g4
pH4(4) = 0; % WNV induced death rate host group 4, gamma4
%pH4(5) = (0.3-log(0.6))/365; % annual per capita birth rate host group 4, Lambda4
if seasons==2
    pH4(5) = 4*(0.3-log(0.6))/365; %4x annual birth rate, assumes all births occur is a three month breeding season.
else
    pH4(5) = 0;
end
pH4(6) = -log(0.6)/365; % natural death rate host group 4, mu_h4
%pH4(7) = 0.0030; % Less Urban carrying capacity host group 4, c_h4
%pH4(8) = 0.0030; % More Urban carrying capacity host group 4, c_h4
pH4(7) = 2; % Less Urban carrying capacity host group 4, c_h4
pH4(8) = 1; % More Urban carrying capacity host group 4, c_h4
pH4(9) = 0.1; % biting preference host group 4, alpha4
%pH4(9) = 1; 
pH4(10) = 1; % prob. horizontal transmission host group 4, p_hh4    %p_hh values all set to 1

%dead-end hosts
pH5 = zeros(1,10);
pH5(1) = 0; % host-to-mosquito transmission probability host group 5, p_hm5
pH5(2) = 0; % maximal horizontal transmission rate host group 5, omega5
pH5(3) = 0; % WNV recovery rate host group 5, g5
pH5(4) = 0; % WNV induced death rate host group 5, gamma5
pH5(5) = 0; % per capita birth rate host group 5, Lambda5
pH5(6) = 0; % natural death rate host group 5, mu_h5
pH5(7) = 70; % Less Urban carrying capacity host group 5, c_h5
pH5(8) = 70; % More Urban carrying capacity host group 5, c_h5
pH5(9) = 0.05; % biting preference host group 5, alpha5
pH5(10) = 0; % prob. horizontal transmission host group 5, p_hh5    %p_hh values all set to 1
 

% --------------------------
% Conditional modifications
% --------------------------

% Disable horizontal transmission for host 1 if transmit == 1
if transmit == 1
    pH1(2) = 0;
end

% Disable host-to-mosquito transmission for host 3 (invulnerable) in case 34
if transmit == 34
    pH3(1) = 0;
end

% Disable host-to-mosquito transmission for host 4 (diluter) in case 34 or 4
if transmit == 34 || transmit == 4
    pH4(1) = 0;
end

% Modify biting preference (alpha) and/or carrying capacity (c_h) under other urban scenarios
switch urban
    case 3
        if tag_input == 1
            pH3(7) = pH3(7) / factor; % carrying capacity
        elseif tag_input == 2
            pH3(9) = pH3(9) / factor; % biting preference
        end
    case 4
        if tag_input == 1
            pH4(7) = pH4(7) * factor;
        elseif tag_input == 2
            pH4(9) = pH4(9) * factor;
        end
    case 5
        pH4(9) = 1;                % Set biting preference of host 4 to 1
        pH1(9) = 1/factor;         % Scale down host 1's biting preference
        pH2(9) = (1/factor) * pH2(9); % Scale down host 2's biting preference
        pH3(9) = (1/factor) * pH3(9); % Scale down host 3's biting preference
end

hostParams = [pH1; pH2; pH3; pH4; pH5];


end