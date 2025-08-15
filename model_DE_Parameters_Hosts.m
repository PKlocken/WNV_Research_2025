function hostParams = model_DE_Parameters_Hosts(urban, transmit, varargin)

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

% --------------------------
% Set default host parameters
% --------------------------

% Host type 1: vulnerable amplifiers with horizontal transmission
pH1 = zeros(1,10);
pH1(1) = 0.68; % host-to-mosquito transmission probability host group 1, p_hm1
pH1(3) = (1/3)*((1/4.7)+(log(0.6)/365))+log(0.6)/365; % WNV recovery rate host group 1, g1
pH1(4) = (1/4.7)+(log(0.6)/365); % WNV induced death rate host group 1, gamma1
pH1(5) = (0.15-log(0.6))/365; % per capita birth rate host group 1, Lambda1
pH1(6) = -log(0.6)/365;       % natural death rate host group 1, mu_h1
pH1(2) = 3.5*(pH1(6)+pH1(4)+pH1(3)); % maximal horizontal transmission rate host group 1, omega1
% parameter comments
pH1(7) = 1;    % carrying capacity (host group 1) in less urban, c_h1
pH1(8) = 0.25; % carrying capacity (host group 1) in more urban, c_h1
pH1(9) = 1;    % biting preference host group 1, alpha1
pH1(10)= 1;    % probability horizontal transmission host group 1, p_hh1

% Host type 2: vulnerable amplifiers without horizontal transmission
pH2 = zeros(1,10);
pH2(1) = 0.53; % p_hm2
pH2(2) = 0;    % omega2
pH2(3) = (1/4.7)+2*(log(0.6)/365); % g2 – WNV recovery rate
pH2(4) = (1/4.7)+(log(0.6)/365);   % gamma2 – induced death rate
pH2(5) = (0.3-log(0.6))/365;       % Lambda2 – birth rate
pH2(6) = -log(0.6)/365;            % mu_h2 – natural death rate
pH2(7) = 2; % c_h2 less urban
pH2(8) = 5; % c_h2 more urban
pH2(9) = 0.1; % alpha2 – biting preference
pH2(10)= 1;   % p_hh2

% Host type 3: invulnerable amplifiers
pH3 = zeros(1,10);
pH3(1) = 0.36; % p_hm3
pH3(2) = 0;    % omega3
pH3(3) = 0.3333; % g3 – recovery rate
pH3(4) = 0;      % gamma3 – no induced death
pH3(5) = (0.1570-log(0.6))/365; % Lambda3
pH3(6) = -log(0.6)/365;         % mu_h3
pH3(7) = 2; % c_h3 less urban
pH3(8) = 1; % c_h3 more urban
pH3(9) = 1; % alpha3
pH3(10)= 1; % p_hh3

% Host type 4: diluters
pH4 = zeros(1,10);
pH4(1) = 0.1;  % p_hm4
pH4(2) = 0;    % omega4
pH4(3) = 1;    % g4 – total recovery
pH4(4) = 0;    % gamma4 – no death
pH4(5) = (0.3-log(0.6))/365; % Lambda4
pH4(6) = -log(0.6)/365;      % mu_h4
pH4(7) = 2; % c_h4 less urban
pH4(8) = 1; % c_h4 more urban
pH4(9) = 0.1; % alpha4
pH4(10)= 1;   % p_hh4

% Dead-end hosts
pH5 = zeros(1,10);
pH5(1) = 0;   % p_hm5
pH5(2) = 0;   % omega5
pH5(3) = 0;   % g5
pH5(4) = 0;   % gamma5
pH5(5) = 0;   % Lambda5
pH5(6) = 0;   % mu_h5
pH5(7) = 70;  % c_h5 less urban
pH5(8) = 70;  % c_h5 more urban
pH5(9) = 0.05; % alpha5
pH5(10)= 0;    % p_hh5


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

% Final output matrix
hostParams = [pH1; pH2; pH3; pH4; pH5];

end
