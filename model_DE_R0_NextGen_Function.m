function [R0, MH_ratio, MH_ratio_a, v_R0] = model_DE_R0_NextGen_Function(N, urban, transmit, extraInput)
% R0 Solver Using Next Gen Matrix Determinant

% Inputs:
%   N               - Number of host types
%   urban           - Urban scenario (1 = less urban, 2 = more urban, 3 = amplify, 4 = dilute, 5=make4 most preferred)
%   transmit        - Transmission toggles (e.g., 4, 34, 1, 0)
%   varargin        - Required if urban > 2 (deppends on urban value e.g., if urban =3,4 then 1 = change density, 2 = change preference)
%
% Outputs:
%   R0              - Basic reproduction number
%   MH_ratio        - Effective mosquito-to-host ratio
%   MH_ratio_a      - Actual mosquito-to-host ratio
%   v_R0            - Eigenvector associated with R0

%-------------------------

if nargin < 4
    extraInput = []; % Default if not provided
end

% Vector Parameters
pV = model_DE_Parameters_Vectors();
rs = pV(1); ri = pV(2); phi = pV(3); qs = pV(4); qi = pV(5);
m_E = pV(6); m_L = pV(7); mu_L = pV(8); mu_V = pV(9);
b = pV(10); c_L = pV(11); k_L = pV(12); p_mh = pV(13);

% Derived Vector Parameters
d_l = ((rs * m_L * qs / mu_V) - mu_L - m_L);
tau = m_L + mu_L + d_l;
M_star = m_L * c_L / mu_V; % DFE susceptible adult mosquitoes

% Vector next-gen matrix blocks
W_M = [m_E 0 0 0; -m_E*qi*phi tau 0 0; 0 0 k_L+mu_V 0; 0 -m_L -k_L mu_V];
F_MM = [0 0 0 ri; 0 0 0 0; 0 0 0 0; 0 0 0 0];

% Retrieve host parameters (with urban and transmission effects)
if urban <= 2
    hostParams = model_DE_Parameters_Hosts(urban, transmit);
else
    hostParams = model_DE_Parameters_Hosts(urban, transmit, extraInput);
end

% Initialize vectors for host matrices
s = [];
omega = [];
f_MH = [];
f_HM = [];

% Totals
NH = 0;
NHa = 0;

% Host contributions to matrices
for j = 1:N
    params = hostParams(j, :);

    p_hm = params(1);
    omega_j = params(2);
    g = params(3);
    gamma = params(4);
    Lambda = params(5);
    mu_h = params(6);

    % Urban/rural carrying capacity
    if urban == 2
        c_h = params(8); % more urban
    else
        c_h = params(7); % less urban or default
    end

    alpha = params(9); % biting preference

    d_h = Lambda - mu_h;
    sj = gamma + g + mu_h + d_h;

    % Add to matrix terms
    s = [s sj];
    omega = [omega omega_j];
    f_MH = [f_MH; b * p_mh * alpha * c_h];
    f_HM = [f_HM b * p_hm * alpha * M_star];

    % Update host totals
    NH = NH + alpha * c_h;
    NHa = NHa + c_h;
end

% Add dead-end host (row 5)
params = hostParams(5, :);
c_h_dead = (urban == 2) * params(8) + (urban ~= 2) * params(7);
alpha_dead = params(9);

NH = NH + alpha_dead * c_h_dead;
NHa = NHa + c_h_dead;

% Ratios
MH_ratio = M_star / NH;   % effective
MH_ratio_a = M_star / NHa; % actual

% Build matrices
W_H = diag(s);
F_HH = diag(omega);
F_HM = zeros(4, N);
F_MH = zeros(N, 4);
F_HM(3, :) = f_HM / NH;
F_MH(:, 4) = f_MH / NH;

% Assemble full matrices
FWinv = [F_HH / W_H, F_MH / W_M;
         F_HM / W_H, F_MM / W_M];

% Compute spectral radius of Next Gen matrix
[eigVecs, eigVals] = eig(FWinv);
eigVals_diag = diag(eigVals);
R0 = max(real(eigVals_diag));

[~, idx] = max(real(eigVals_diag));
v_R0 = eigVecs(:, idx);
v_R0 = v_R0 / max(abs(v_R0));  % Normalize

% Display (optional)
display(R0);
display(MH_ratio);
display(MH_ratio_a);
display(v_R0);
end
