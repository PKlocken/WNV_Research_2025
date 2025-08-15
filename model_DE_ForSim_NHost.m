function [dxdt] = model_DE_ForSim_NHost(t,x,N,urban, transmit, varargin)

% Inputs:
%   N              - Number of host types
%   urban           - Urban scenario (1 = less urban, 2 = more urban, 3 = amplify, 4 = dilute, 5=make4 most preferred)
%   transmit        - Transmission toggles (e.g., 4, 34, 1, 0)
%   varargin        - Required if urban > 2 (deppends on urban value e.g., if urban =3,4 then 1 = change density, 2 = change preference)

% Outputs:
%   dxdt           - the derivatives of the continuous state variables.

%-------------------------

if urban > 2 && urban < 5
    if isempty(varargin)
        error('Additional input required for urban > 2 cases.');
    end
    tag_input = varargin{1};
    factor = 100;
end

% Vector Parameter Data
pV = model_DE_Parameters_Vectors();

% Vector Parameters for Model
rs = pV(1);
ri = pV(2);
phi = pV(3);
qs = pV(4);
qi = pV(5);
m_E = pV(6);
m_L = pV(7);
mu_L = pV(8);
mu_V = pV(9);
b = pV(10);
c_L = pV(11);
kl = pV(12);
p_mh = pV(13);
d_l = ((rs*m_L*qs/mu_V)-mu_L-m_L); % density-dependent larval death

% Vector State variables
Es = x(1);
Ei = x(2);
Ls = x(3);
Li = x(4);
Vs = x(5);
Ve = x(6);
Vi = x(7);

%-------------------------

% Control Parameter Data
pU = model_DE_Parameters_Control(1, 10);
km1 = pU(1);
km2 = pU(4);

% Control State Variables
Ul = 0;
Ua = 0;

%-------------------------

% Host Parameter Data — NOW receives urban, transmit, and varargin
hostParams = model_DE_Parameters_Hosts(urban, transmit, varargin{:});

% Initialize Host Compartments
NH = zeros(N,1);
Y = zeros(N,1);
alpha = zeros(N,1);

dH = [];
dVs = m_L * Ls;
dVe = 0;

% Loop through host types
for j = 1:N
    pH = hostParams(j, :);

    alpha(j) = pH(9); % biting preference

    Hs = x(8 + 3*(j-1));
    Hi = x(9 + 3*(j-1));
    Hr = x(10 + 3*(j-1));
    NH(j) = Hs + Hi + Hr;

    Y(j) = alpha(j) * NH(j);
end

NY = sum(Y);

% Dead-end host group (5)
params = hostParams(5, :);
if urban == 2
    c_h = params(8); % More urban
else
    c_h = params(7); % Less urban
end
alphaD = params(9);
NY = NY + alphaD * c_h;

for j = 1:N
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

    NHj = NH(j);
    Hs = x(8 + 3*(j-1));
    Hi = x(9 + 3*(j-1));
    Hr = x(10 + 3*(j-1));

    % Host ODEs
    dHs = Lambda * NHj - b * p_mh * Vi * alpha(j) * Hs / NY - omega * p_hh * Hi * Hs / NHj - d_h * NHj * Hs / c_h - mu_h * Hs;
    dHi = b * p_mh * Vi * alpha(j) * Hs / NY + omega * p_hh * Hi * Hs / NHj - (gamma + g) * Hi - d_h * NHj * Hi / c_h - mu_h * Hi;
    dHr = g * Hi - d_h * NHj * Hr / c_h - mu_h * Hr;

    dH = [dH; dHs; dHi; dHr];

    % Vector updates
    dVs = dVs - b * p_hm * Vs * alpha(j) * Hi / NY;
    dVe = dVe + b * p_hm * Vs * alpha(j) * Hi / NY;
end

% Final vector equations
dVs = dVs - mu_V * Vs - km2 * Vs * Ua;
dVe = dVe - kl * Ve - mu_V * Ve - km2 * Ve * Ua;

% Other vector dynamics
dEs = rs * (Vs + Ve) - m_E * Es;
dEi = ri * Vi - m_E * Ei;
dLs = m_E * qs * Es + m_E * qi * (1 - phi) * Ei - mu_L * Ls - m_L * Ls - d_l * Ls * (Ls + Li) / c_L - km1 * Ls * Ul;
dLi = m_E * qi * phi * Ei - mu_L * Li - m_L * Li - d_l * Li * (Li + Ls) / c_L - km1 * Li * Ul;
dVi = m_L * Li + kl * Ve - mu_V * Vi - km2 * Vi * Ua;

% Combine derivatives
dxdt = [dEs; dEi; dLs; dLi; dVs; dVe; dVi; dH];

end
