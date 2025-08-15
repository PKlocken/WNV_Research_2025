function [dxdt] = model_DE_ForSim_NHost_seasons(t,x,N,urban,transmit, seasons, varargin)

% N is number of host types
% urban = 1 means less urban
% urban = 2 means more urban
% This file gives the derivatives of the continuous state variables.

%-------------------------

% Vector Parameter Data
pV = model_DE_Parameters_Vectors_seasons(seasons);

% Vector Parameters for Model
%assume mosquitoes breed for half the year April-September
rs = pV(1); % egg laying rate of S and E mosquitoes
ri = pV(2); % egg laying rate of I mosquitoes
phi = pV(3); % fraction of eggs born to infected mothers that are infected (fraction of eggs infected)
qs = pV(4); % fraction of eggs from uninfected mosquitoes that hatch
qi = pV(5); % fraction of eggs laid to infected mosquitoes that hatch
m_E = pV(6); % hatch rate 
m_L = pV(7); % larval maturation rate
mu_L = pV(8); % larval death rate
mu_V = pV(9); % adult death rate
b = pV(10); % mosquito biting rate
c_L = pV(11); % mosquito carrying capacity (larval)
kl = pV(12); % disease progression in mosquitoes (1/latency period)
p_mh = pV(13); % mosquito-to-host transmission
d_l=((rs*m_L*qs/mu_V)-mu_L-m_L); % density-dependent death rate for larvae

% Vector State variables
Es = x(1); % eggs laid by susceptible and exposed mothers
Ei = x(2); % eggs laid by infected mothers
Ls = x(3); % susceptible larvae
Li = x(4); % infected larvae
Vs = x(5); % susceptible vectors
Ve = x(6); % exposed vectors
Vi = x(7); % infected vectors

%-------------------------

% Control Parameter Data
pU = model_DE_Parameters_Control(1, 10);

% Control Parameters for Model
km1 = pU(1); % max rate at which larvicide kills larvae
km2 = pU(4); % max rate at which adulticide kills adult vectors

% Control State Variables
Ul = 0; % larvicide
Ua = 0; % adulticide

%-------------------------

% Host Parameter Data
hostParams = model_DE_Parameters_Hosts_seasons(seasons, urban, transmit, varargin{:});

% Initialize Host Compartments
NH = zeros(N,1); % jth component gives number of hosts of type j
Y = zeros(N,1); % jth component gives scaled number of hosts of type j
alpha = zeros(N,1); % jth component gives alpha for hosts of type j

dH = [];
dVs = m_L*Ls;
dVe = 0;

% Loop through host types dynamically
for j = 1:N
    % Retrieve host-specific parameters
    pH = hostParams(j, :);
    alpha(j) = pH(9); % biting preference

    % State Variables for current host type
    Hs = x(8 + 3*(j-1)); % susceptible hosts
    Hi = x(9 + 3*(j-1)); % infected hosts
    Hr = x(10 + 3*(j-1)); % recovered hosts
    NH(j) = Hs + Hi + Hr; % total hosts of type j

    % Update totals with biting preference weight
    Y(j) = alpha(j) * NH(j);
end

NY=sum(Y);

%dead-end hosts
params = hostParams(5, :);
if urban == 2
        c_h = params(8); % More urban carrying capacity
    else
        c_h = params(7); % Less urban carrying capacity
end
alphaD = params(9); % Biting preference
NY=NY+alphaD*c_h;

for j=1:N

    % Retrieve host-specific parameters
    pH = hostParams(j, :);

    p_hm = pH(1);
    omega = pH(2);
    g = pH(3);
    gamma = pH(4);
    Lambda = pH(5);
    mu_h = pH(6);

    switch urban
        case 2
            c_h = pH(8);
        otherwise
            c_h = pH(7);
    end

    p_hh = pH(10);
    d_h = Lambda - mu_h;

    NHj=NH(j);

    % State Variables for current host type
    Hs = x(8 + 3*(j-1)); % susceptible hosts
    Hi = x(9 + 3*(j-1)); % infected hosts
    Hr = x(10 + 3*(j-1)); % recovered hosts
    
    % ODEs for hosts of type j
    dHs = Lambda * NHj - b * p_mh * Vi * alpha(j) * Hs / NY - omega * p_hh * Hi * Hs / NHj - d_h * NHj * Hs / c_h - mu_h * Hs;
    dHi = b * p_mh * Vi * alpha(j) * Hs / NY + omega * p_hh * Hi * Hs / NHj - (gamma + g) * Hi - d_h * NHj * Hi / c_h - mu_h * Hi;
    dHr = g * Hi - d_h * NHj * Hr / c_h - mu_h * Hr;

    % Append derivatives for current host type
    dH = [dH; dHs; dHi; dHr];

    % Update vector ODEs
    dVs = dVs - b * p_hm * Vs * alpha(j) * Hi / NY;
    dVe = dVe + b * p_hm * Vs * alpha(j) * Hi / NY;
end

% Finalize vector ODEs

dVs = dVs - mu_V * Vs - km2 * Vs * Ua;
dVe = dVe - kl * Ve - mu_V * Ve - km2 * Ve * Ua;

% Other Vector ODEs
dEs = rs * (Vs + Ve) - m_E * Es;
dEi = ri * Vi - m_E * Ei;
dLs = m_E * qs * Es + m_E * qi * (1 - phi) * Ei - mu_L * Ls - m_L * Ls - d_l * Ls * (Ls + Li) / c_L - km1 * Ls * Ul;
dLi = m_E * qi * phi * Ei - mu_L * Li - m_L * Li - d_l * Li * (Li + Ls) / c_L - km1 * Li * Ul;
dVi = m_L * Li + kl * Ve - mu_V * Vi - km2 * Vi * Ua;

% Compile derivatives
dxdt = [dEs; dEi; dLs; dLi; dVs; dVe; dVi; dH];

end
