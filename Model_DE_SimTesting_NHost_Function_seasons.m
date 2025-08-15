function [hostData_long_final, Max_Prevalence, R0, MH_ratio, MH_ratio_a, v_R0] = Model_DE_SimTesting_NHost_Function_seasons(N, urban, transmit, infection_origin, num_years, savePath, extraInput)

% N = 4; % Number of host types
% urban = 1; % Urban parameter (1 = less urban, 2 = more urban)
% % 
% % 
% % transmit = 4; %P_{HM}=0 for host 4
% % transmit = 34; %P_{HM}=0 for hosts 3 and 4
% % transmit = 1; % with omega=0 for host 1
% transmit = 0; % All transmissions active


% If savePath is not provided or is empty, create it with timestamp
if nargin < 6 || isempty(savePath)
    timestamp = string(datetime('now', 'Format', 'yyyy_MM_dd_HHmmss'));
    savePath = fullfile('Plots', "Run_" + timestamp);
    mkdir(savePath);
elseif ~exist(savePath, 'dir')
    % savePath provided but folder does not exist → create it
    mkdir(savePath);
end

if urban > 2 && urban < 5
    if isempty(extraInput)
        error('extraInput is required when urban = 3 or 4.');
    end
    tag_input = extraInput;
    switch tag_input
        case 1
            tag = 'cd'; % capacity dilution
        case 2
            tag = 'cb'; % biting preference
        otherwise
            tag = 'cn'; % control/no-change
    end
else
    tag_input = [];
    tag = 'cn'; % default tag for urban <= 2 or urban >= 5
    extraInput = [];
end

T_short_initial = 0;
T_short_final = 250;
T_long_initial = 250;
T_long_final = num_years*365;
tspan = [T_short_initial, T_long_final];

% Vector Parameter Data
pV = model_DE_Parameters_Vectors();

% Vector Parameters for Model
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

% Vector Disease-Free Equilibrium
Es_DFE = c_L*m_L*rs/(mu_V*m_E);
Ls_DFE = c_L;
Vs_DFE = c_L*m_L/mu_V;
origin_label = [];
host_originator = 0;

switch infection_origin % Initial conditions for vectors (Es, Ei, Ls, Li, Vs, Ve, Vi)
    case 'e' % infection originates with eggs
        x0 = [0.99*Es_DFE; 0.01*Es_DFE; Ls_DFE; 0; Vs_DFE; 0; 0];
        origin_label = "Eggs";
    case 'l' % infection originates with larvea
        x0 = [Es_DFE; 0; 0.99*Ls_DFE; 0.01*Ls_DFE; Vs_DFE; 0; 0];
        origin_label = "Larvea";
    case 'M' % infection originates with adult exposed mosquito
        x0 = [Es_DFE; 0; Ls_DFE; 0; 0.99*Vs_DFE; 0.01*Vs_DFE; 0];
        origin_label = "Exposed Mosquitoes";
    case 'm' % infection originates with adult infected mosquito
        %x0 = [Es_DFE; 0; Ls_DFE; 0; 0.99*Vs_DFE; 0; 0.01*Vs_DFE];
        x0 = [Es_DFE; 0; Ls_DFE; 0; 0.999*Vs_DFE; 0; 0.001*Vs_DFE];
        origin_label = "Infected Mosquitoes";
        %infection_origin = 'm';
    otherwise
        x0 = [Es_DFE; 0; Ls_DFE; 0; Vs_DFE; 0; 0];
        host_originator = str2double(infection_origin);
end

% Load host parameters with urban/transmission logic centralized
if urban <= 2
    hostParams = model_DE_Parameters_Hosts(urban, transmit);
else
    hostParams = model_DE_Parameters_Hosts(urban, transmit, extraInput);
end

% Add initial conditions for hosts dynamically
for j = 1:N
    % Retrieve host-specific parameters

    pH = hostParams(j, :);
    if urban == 2
        c_h = pH(8);
    else
        c_h = pH(7);
    end

    Hs_DFE = c_h; % Disease-Free Equilibrium, host type j

    if host_originator == j
        x0 = [x0; 0.99*Hs_DFE; 0.01*Hs_DFE; 0]; % Initial conditions for each host type (Hs, Hi, Hr)
    else
        x0 = [x0; Hs_DFE; 0; 0]; % Initial conditions for each host type (Hs, Hi, Hr)
    end
end


t=[];
x=[];
x01=x0;


if urban <= 2
    % Calculate Reproduction Number
    [R0, MH_ratio, MH_ratio_a, v_R0] = model_DE_R0_NextGen_Function(N, urban, transmit, []);
    for i=1:num_years
        tspan1 = (i-1)*365+[0 365/12]; % only mosquuitoes breed
        tspan2 = (i-1)*365+[365/12 4*365/12]; %birds and mosquitoes breed
        tspan3 = (i-1)*365+[4*365/12 7*365/12]; %only mosquitoes breed
        tspan4 = (i-1)*365+[7*365/12 365];%no breeding

        % Solve ODE System
        [t1, x1] = ode45(@(t, x) model_DE_ForSim_NHost_seasons(t, x, N, urban, transmit,1), tspan1, x01);
        x02=x1(length(t1),:);
        [t2, x2] = ode45(@(t, x) model_DE_ForSim_NHost_seasons(t, x, N, urban, transmit,2), tspan2, x02);
        x03=x2(length(t2),:);
        [t3, x3] = ode45(@(t, x) model_DE_ForSim_NHost_seasons(t, x, N, urban, transmit,3), tspan3, x03);
        x04=x3(length(t3),:);
        [t4, x4] = ode45(@(t, x) model_DE_ForSim_NHost_seasons(t, x, N, urban, transmit,4), tspan4, x04);
        x01=x4(length(t4),:);
        t=[t; t1; t2; t3; t4];
        x=[x; x1; x2; x3; x4];
    end
else
    % Calculate Reproduction Number
    [R0, MH_ratio, MH_ratio_a, v_R0] = model_DE_R0_NextGen_Function(N, urban, transmit, extraInput);
    for i=1:num_years
        tspan1 = (i-1)*365+[0 365/12]; % only mosquuitoes breed
        tspan2 = (i-1)*365+[365/12 4*365/12]; %birds and mosquitoes breed
        tspan3 = (i-1)*365+[4*365/12 7*365/12]; %only mosquitoes breed
        tspan4 = (i-1)*365+[7*365/12 365];%no breeding

        % Solve ODE System
        [t1, x1] = ode45(@(t, x) model_DE_ForSim_NHost_seasons(t, x, N, urban, transmit,1,extraInput), tspan1, x01);
        x02=x1(length(t1),:);
        [t2, x2] = ode45(@(t, x) model_DE_ForSim_NHost_seasons(t, x, N, urban, transmit,2,extraInput), tspan2, x02);
        x03=x2(length(t2),:);
        [t3, x3] = ode45(@(t, x) model_DE_ForSim_NHost_seasons(t, x, N, urban, transmit,3,extraInput), tspan3, x03);
        x04=x3(length(t3),:);
        [t4, x4] = ode45(@(t, x) model_DE_ForSim_NHost_seasons(t, x, N, urban, transmit,4,extraInput), tspan4, x04);
        x01=x4(length(t4),:);
        t=[t; t1; t2; t3; t4];
        x=[x; x1; x2; x3; x4];
    end
end

% Extract time ranges
short_time_range = (t >= T_short_initial & t <= T_short_final); % Short time range
long_time_range = (t >= T_long_initial & t <= T_long_final); % Long time range

% Extract data for short time range
t_short = t(short_time_range);
x_short = x(short_time_range, :);

% Extract data for long time range
t_long = t(long_time_range);
x_long = x(long_time_range, :);

% Extract vector data, short time range
Es_short = x_short(:, 1);
Ei_short = x_short(:, 2);
Ls_short = x_short(:, 3);
Li_short = x_short(:, 4);
Vs_short = x_short(:, 5);
Ve_short = x_short(:, 6);
Vi_short = x_short(:, 7);
Mosquito_Prevalence = Vi_short./(Vi_short+Ve_short+Vs_short);
Max_Prevalence = max(Mosquito_Prevalence);

% Extract vector data, long time range
Es_long = x_long(:, 1);
Ei_long = x_long(:, 2);
Ls_long = x_long(:, 3);
Li_long = x_long(:, 4);
Vs_long = x_long(:, 5);
Ve_long = x_long(:, 6);
Vi_long = x_long(:, 7);


% Extract host data for both ranges
hostData_short = x_short(:, 8:end);
hostData_long = x_long(:, 8:end);

%[r,c] = size(x_long);
hostData_long_final = x_long(length(t_long), 8:end);
hostData_short_final = x_short(length(t_short), 8:end);
% length(long_time_range)
% length(t_long)
% size(x_long)
display(Max_Prevalence);
% return

total_short_final = sum(reshape(hostData_short_final, 3, []), 1);
total_long_final  = sum(reshape(hostData_long_final, 3, []), 1);

disp('Total final densities (short time):');
disp(total_short_final);

disp('Total final densities (long time):');
disp(total_long_final);



%================================================================================================================

figureTitleFontSize = 25;
plotTitleFontSize = 25;
plotIndexFontSize = 20;

% Combined Infected Compartments Plots, Short Time
figure('Position', get(0, 'Screensize'));
tlo = tiledlayout(2, 1); % Default padding and spacing

% Plot all infected hosts
ax1 = nexttile;
hold(ax1, 'on');
for j = 1:N
    Hi = hostData_short(:, 3*(j-1) + 2);
    plot(t_short, Hi, 'LineWidth', 4);
end
title(ax1, 'Infected Hosts', 'FontSize', plotTitleFontSize);
xlabel(ax1, 'Time (days)', 'FontSize', plotIndexFontSize);
ylabel(ax1, 'Population (ind/ha)', 'FontSize', plotIndexFontSize);
legend(ax1, arrayfun(@(j) ['Type ' num2str(j)], 1:N, 'UniformOutput', false), 'Location', 'best', 'FontSize', plotIndexFontSize);
grid(ax1, 'on');

% Plot infected and exposed mosquitoes
ax2 = nexttile;
hold(ax2, 'on');
plot(t_short, Ve_short, 'm', 'LineWidth', 4);
plot(t_short, Vi_short, 'r', 'LineWidth', 4);
title(ax2, 'Infected Vectors', 'FontSize', plotTitleFontSize);
xlabel(ax2, 'Time (days)', 'FontSize', plotIndexFontSize);
ylabel(ax2, 'Population (ind/ha)', 'FontSize', plotIndexFontSize);
legend(ax2, {'Exposed Mosquitoes', 'Infected Mosquitoes'}, 'Location', 'best', 'FontSize', plotIndexFontSize);
grid(ax2, 'on');

% if urban == 2
%     switch transmit
%         case 4
%             sgtitle(sprintf('Infected Compartments, More Urban, No Host-to-Mosquito Transmission for Type 4, R0=%.2f, Origin: %s, Seasonal', R0, origin_label), 'FontSize', figureTitleFontSize);
%         case 34
%             sgtitle(sprintf('Infected Compartments, More Urban, No Host-to-Mosquito Transmission for Type 3 or 4, R0=%.2f, Origin: %s, Seasonal', R0, origin_label), 'FontSize', figureTitleFontSize);
%         case 1
%             sgtitle(sprintf('Infected Compartments, More Urban, No Host-to-Host Transmission for Type 1, R0=%.2f, Origin: %s, Seasonal', R0, origin_label), 'FontSize', figureTitleFontSize);
%         otherwise
%             sgtitle(sprintf('Infected Compartments, More Urban, All Transmissions, R0=%.2f, Origin: %s, Seasonal', R0, origin_label), 'FontSize', figureTitleFontSize);
%     end
% else 
%     switch transmit
%         case 4
%             sgtitle(sprintf('Infected Compartments, Less Urban, No Host-to-Mosquito Transmission for Type 4, R0=%.2f, Origin: %s, Seasonal', R0, origin_label), 'FontSize', figureTitleFontSize);
%         case 34
%             sgtitle(sprintf('Infected Compartments, Less Urban, No Host-to-Mosquito Transmission for Type 3 or 4, R0=%.2f, Origin: %s, Seasonal', R0, origin_label), 'FontSize', figureTitleFontSize);
%         case 1
%             sgtitle(sprintf('Infected Compartments, Less Urban, No Host-to-Host Transmission for Type 1, R0=%.2f, Origin: %s, Seasonal', R0, origin_label), 'FontSize', figureTitleFontSize);
%         otherwise
%             sgtitle(sprintf('Infected Compartments, Less Urban, All Transmissions, R0=%.2f, Origin: %s, Seasonal', R0, origin_label), 'FontSize', figureTitleFontSize);
%     end
% end


% Build figure caption text (replace sgtitle)
% if urban == 2
%     switch transmit
%         case 4
%             figCaption = sprintf('Infected Compartments, More Urban, No Host-to-Mosquito Transmission for Type 4, R0=%.2f, Origin: %s, Seasonal', R0, origin_label);
%         case 34
%             figCaption = sprintf('Infected Compartments, More Urban, No Host-to-Mosquito Transmission for Type 3 or 4, R0=%.2f, Origin: %s, Seasonal', R0, origin_label);
%         case 1
%             figCaption = sprintf('Infected Compartments, More Urban, No Host-to-Host Transmission for Type 1, R0=%.2f, Origin: %s, Seasonal', R0, origin_label);
%         otherwise
%             figCaption = sprintf('Infected Compartments, More Urban, All Transmissions, R0=%.2f, Origin: %s, Seasonal', R0, origin_label);
%     end
% else
%     switch transmit
%         case 4
%             figCaption = sprintf('Infected Compartments, Less Urban, No Host-to-Mosquito Transmission for Type 4, R0=%.2f, Origin: %s, Seasonal', R0, origin_label);
%         case 34
%             figCaption = sprintf('Infected Compartments, Less Urban, No Host-to-Mosquito Transmission for Type 3 or 4, R0=%.2f, Origin: %s, Seasonal', R0, origin_label);
%         case 1
%             figCaption = sprintf('Infected Compartments, Less Urban, No Host-to-Host Transmission for Type 1, R0=%.2f, Origin: %s, Seasonal', R0, origin_label);
%         otherwise
%             figCaption = sprintf('Infected Compartments, Less Urban, All Transmissions, R0=%.2f, Origin: %s, Seasonal', R0, origin_label);
%     end
% end

% Add caption below all subplots using tiledlayout's xlabel
%xlabel(tlo, figCaption, 'FontSize', figureTitleFontSize, 'Interpreter', 'none', 'HorizontalAlignment', 'center');

% Save figure as before
switch transmit
    case 4
        file_base = sprintf('infected_compartments_T0=%.0f_Tf=%.0f_N=%.0f_urban=%.0f_no4_seasons__R0=%.2f_%s', T_short_initial, T_short_final, N, urban, R0, infection_origin);
    case 34 
        file_base = sprintf('infected_compartments_T0=%.0f_Tf=%.0f_N=%.0f_urban=%.0f_no34_seasons__R0=%.2f_%s', T_short_initial, T_short_final, N, urban, R0, infection_origin);
    case 1
        file_base = sprintf('infected_compartments_T0=%.0f_Tf=%.0f_N=%.0f_urban=%.0f_no1_seasons__R0=%.2f_%s', T_short_initial, T_short_final, N, urban, R0, infection_origin);
    otherwise
        file_base = sprintf('infected_compartments_T0=%.0f_Tf=%.0f_N=%.0f_urban=%.0f_all_seasons__R0=%.2f_%s', T_short_initial, T_short_final, N, urban, R0, infection_origin);
end

file_path_eps = fullfile(savePath, [file_base, '.eps']);
file_path_fig = fullfile(savePath, [file_base, '.fig']);

saveas(gcf, file_path_eps, 'epsc');
saveas(gcf, file_path_fig, 'fig');



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


% Combined Infected Compartments Plot, Long Time
figure('Position', get(0, 'Screensize'));
tiledlayout(2, 1, 'Padding', 'none', 'TileSpacing', 'none'); % Compact layout

% Plot all infected hosts
subplot(2, 1, 1);
hold on;
for j = 1:N
    Hi = hostData_long(:, 3*(j-1) + 2);
    plot(t_long, Hi, 'LineWidth', 4); % Plot each infected group
end
xlim([T_long_initial, max(t_long)]);
ylim tight
title('Infected Hosts', 'FontSize', plotTitleFontSize);
xlabel('Time (days)', 'FontSize', plotIndexFontSize);
ylabel('Population (ind/ha)', 'FontSize', plotIndexFontSize);
legend(arrayfun(@(j) ['Type ' num2str(j)], 1:N, 'UniformOutput', false), 'Location', 'best', 'FontSize', plotIndexFontSize);
grid on;

% Plot exposed mosquitoes and infected mosquitoes
subplot(2, 1, 2);
plot(t_long, Ve_long, 'm', 'LineWidth', 4); hold on;
plot(t_long, Vi_long, 'r', 'LineWidth', 4);
xlim([T_long_initial, max(t_long)]);
ylim tight
title('Infected Vectors', 'FontSize', plotTitleFontSize);
xlabel('Time (days)', 'FontSize', plotIndexFontSize);
ylabel('Population (ind/ha)', 'FontSize', plotIndexFontSize);
legend('Exposed Mosquitoes', 'Infected Mosquitoes', 'Location', 'best', 'FontSize', plotIndexFontSize);
grid on;

%Figure Title
% if urban == 2
%     switch transmit
%         case 4
%             sgtitle(sprintf('Infected Compartments, More Urban, No Host-to-Mosquito Transmission for Type 4, R0=%.2f, Origin: %s, Seasonal', R0, origin_label), 'FontSize', figureTitleFontSize);
%         case 34
%             sgtitle(sprintf('Infected Compartments, More Urban, No Host-to-Mosquito Transmission for Type 3 or 4, R0=%.2f, Origin: %s, Seasonal', R0, origin_label), 'FontSize', figureTitleFontSize);
%         case 1
%             sgtitle(sprintf('Infected Compartments, More Urban, No Host-to-Host Transmission for Type 1, R0=%.2f, Origin: %s, Seasonal', R0, origin_label), 'FontSize', figureTitleFontSize);
%         otherwise
%             sgtitle(sprintf('Infected Compartments, More Urban, All Transmissions, R0=%.2f, Origin: %s, Seasonal', R0, origin_label), 'FontSize', figureTitleFontSize);
%     end
% else 
%     switch transmit
%         case 4
%             sgtitle(sprintf('Infected Compartments, Less Urban, No Host-to-Mosquito Transmission for Type 4, R0=%.2f, Origin: %s, Seasonal', R0, origin_label), 'FontSize', figureTitleFontSize);
%         case 34
%             sgtitle(sprintf('Infected Compartments, Less Urban, No Host-to-Mosquito Transmission for Type 3 or 4, R0=%.2f, Origin: %s, Seasonal', R0, origin_label), 'FontSize', figureTitleFontSize);
%         case 1
%             sgtitle(sprintf('Infected Compartments, Less Urban, No Host-to-Host Transmission for Type 1, R0=%.2f, Origin: %s, Seasonal', R0, origin_label), 'FontSize', figureTitleFontSize);
%         otherwise
%             sgtitle(sprintf('Infected Compartments, Less Urban, All Transmissions, R0=%.2f, Origin: %s, Seasonal', R0, origin_label), 'FontSize', figureTitleFontSize);
%     end
% end

% Save the figure as a PDF
switch transmit
    case 4
        file_base = sprintf('infected_compartments_T0=%.0f_Tf=%.0f_N=%.0f_urban=%.0f_no4_seasons__R0=%.2f_%s', T_long_initial, T_long_final, N, urban, R0, infection_origin);
    case 34 
        file_base = sprintf('infected_compartments_T0=%.0f_Tf=%.0f_N=%.0f_urban=%.0f_no34_seasons__R0=%.2f_%s', T_long_initial, T_long_final, N, urban, R0, infection_origin);
    case 1
        file_base = sprintf('infected_compartments_T0=%.0f_Tf=%.0f_N=%.0f_urban=%.0f_no1_seasons__R0=%.2f_%s', T_long_initial, T_long_final, N, urban, R0, infection_origin);
    otherwise
        file_base = sprintf('infected_compartments_T0=%.0f_Tf=%.0f_N=%.0f_urban=%.0f_all_seasons__R0=%.2f_%s', T_long_initial, T_long_final, N, urban, R0, infection_origin);
end

% Construct full file paths
file_path_eps = fullfile(savePath, [file_base, '.eps']);
file_path_fig = fullfile(savePath, [file_base, '.fig']);

% Save figures
saveas(gcf, file_path_eps, 'epsc');
saveas(gcf, file_path_fig, 'fig');

%==============================================
%==============================================

% Combined Susceptible/Recovered Host Compartments ShortTime Plots

figure('Position', get(0, 'Screensize'));
tiledlayout(2, 1); % Default padding and spacing

% Plot all susceptible hosts, short time
ax1 = nexttile;
hold(ax1, 'on');
for j = 1:N
    Hs = hostData_short(:, 3*(j-1) + 1);
    plot(t_short, Hs, 'LineWidth', 4); % Plot each susceptible group
end
title('Susceptible Hosts', 'FontSize', plotTitleFontSize);
xlabel('Time (days)', 'FontSize', plotIndexFontSize);
ylabel('Population (ind/ha)', 'FontSize', plotIndexFontSize);
legend(arrayfun(@(j) ['Type ' num2str(j)], 1:N, 'UniformOutput', false), 'Location', 'best', 'FontSize', plotIndexFontSize);
grid(ax1, 'on');

% Plot all recoverd hosts, short time
ax2 = nexttile;
hold(ax2, 'on');
for j = 1:N
    Hr = hostData_short(:, 3*(j-1) + 3);
    plot(t_short, Hr, 'LineWidth', 4); % Plot each recovered group
end
title('Recovered Hosts', 'FontSize', plotTitleFontSize);
xlabel('Time (days)', 'FontSize', plotIndexFontSize);
ylabel('Population (ind/ha)', 'FontSize', plotIndexFontSize);
legend(arrayfun(@(j) ['Type ' num2str(j)], 1:N, 'UniformOutput', false), 'Location', 'best', 'FontSize', plotIndexFontSize);
grid(ax2, 'on');

% Save figure as before
switch transmit
    case 4
        file_base = sprintf('recovered_compartments_T0=%.0f_Tf=%.0f_N=%.0f_urban=%.0f_no4_seasons__R0=%.2f_%s', T_short_initial, T_short_final, N, urban, R0, infection_origin);
    case 34 
        file_base = sprintf('recovered_compartments_T0=%.0f_Tf=%.0f_N=%.0f_urban=%.0f_no34_seasons__R0=%.2f_%s', T_short_initial, T_short_final, N, urban, R0, infection_origin);
    case 1
        file_base = sprintf('recovered_compartments_T0=%.0f_Tf=%.0f_N=%.0f_urban=%.0f_no1_seasons__R0=%.2f_%s', T_short_initial, T_short_final, N, urban, R0, infection_origin);
    otherwise
        file_base = sprintf('recovered_compartments_T0=%.0f_Tf=%.0f_N=%.0f_urban=%.0f_all_seasons__R0=%.2f_%s', T_short_initial, T_short_final, N, urban, R0, infection_origin);
end

file_path_eps = fullfile(savePath, [file_base, '.eps']);
file_path_fig = fullfile(savePath, [file_base, '.fig']);

saveas(gcf, file_path_eps, 'epsc');
saveas(gcf, file_path_fig, 'fig');


%=======================================

% Combined Susceptible/Recovered Host Compartments Long_time Plots

figure('Position', get(0, 'Screensize'));
tiledlayout(2, 1); % Default padding and spacing

% Plot all susceptible hosts, long time
ax1 = nexttile;
hold(ax1, 'on');
for j = 1:N
    Hs = hostData_long(:, 3*(j-1) + 1);
    plot(t_long, Hs, 'LineWidth', 4); % Plot each susceptible group
end
xlim([T_long_initial, max(t_long)]);
ylim tight
title('Susceptible Hosts', 'FontSize', plotTitleFontSize);
xlabel('Time (days)', 'FontSize', plotIndexFontSize);
ylabel('Population (ind/ha)', 'FontSize', plotIndexFontSize);
legend(arrayfun(@(j) ['Type ' num2str(j)], 1:N, 'UniformOutput', false), 'Location', 'best', 'FontSize', plotIndexFontSize);
grid(ax1, 'on');

% Plot all recovered hosts, long time
ax2 = nexttile;
hold(ax2, 'on');
for j = 1:N
    Hr = hostData_long(:, 3*(j-1) + 3);
    plot(t_long, Hr, 'LineWidth', 4); % Plot each recovered group
end
xlim([T_long_initial, max(t_long)]);
ylim tight
title('Recovered Hosts', 'FontSize', plotTitleFontSize);
xlabel('Time (days)', 'FontSize', plotIndexFontSize);
ylabel('Population (ind/ha)', 'FontSize', plotIndexFontSize);
legend(arrayfun(@(j) ['Type ' num2str(j)], 1:N, 'UniformOutput', false), 'Location', 'best', 'FontSize', plotIndexFontSize);
grid(ax2, 'on');

% Save figure as before
switch transmit
    case 4
        file_base = sprintf('noninfected_compartments_T0=%.0f_Tf=%.0f_N=%.0f_urban=%.0f_no4_seasons__R0=%.2f_%s', T_long_initial, T_long_final, N, urban, R0, infection_origin);
    case 34 
        file_base = sprintf('noninfected_compartments_T0=%.0f_Tf=%.0f_N=%.0f_urban=%.0f_no34_seasons__R0=%.2f_%s', T_long_initial, T_long_final, N, urban, R0, infection_origin);
    case 1
        file_base = sprintf('noninfected_compartments_T0=%.0f_Tf=%.0f_N=%.0f_urban=%.0f_no1_seasons__R0=%.2f_%s', T_long_initial, T_long_final, N, urban, R0, infection_origin);
    otherwise
        file_base = sprintf('noninfected_compartments_T0=%.0f_Tf=%.0f_N=%.0f_urban=%.0f_all_seasons__R0=%.2f_%s', T_long_initial, T_long_final, N, urban, R0, infection_origin);
end

file_path_eps = fullfile(savePath, [file_base, '.eps']);
file_path_fig = fullfile(savePath, [file_base, '.fig']);

saveas(gcf, file_path_eps, 'epsc');
saveas(gcf, file_path_fig, 'fig');

% Create a summary .txt file with function outputs


if isnumeric(extraInput)
    extraStr = num2str(extraInput);
else
    extraStr = extraInput;
end

if isempty(extraStr)
    filename = sprintf('SimResults_seasonal_Urban%d_Trans%.2f_%s.txt', urban, transmit, tag);
else
    filename = sprintf('SimResults_seasonal_Urban%d_Trans%.2f_%s_Extra%s.txt', urban, transmit, tag, extraStr);
end

txtFileName = fullfile(savePath, filename);
fid = fopen(txtFileName, 'w');

if fid == -1
    warning('Unable to create summary text file.');
else
    fprintf(fid, 'Simulation Summary\n');
    fprintf(fid, '==================\n\n');
    fprintf(fid, 'R0: %.6f\n', R0);
    fprintf(fid, 'Corresponding Eigenvector: ');
    fprintf(fid, '%.6f ', v_R0);
    fprintf(fid, '\n');
    fprintf(fid, 'Max Mosquito Prevalence: %.6f\n', Max_Prevalence);
    fprintf(fid, 'Mosquito-to-Host Ratio (MH_ratio): %.6f\n', MH_ratio);
    fprintf(fid, 'Adjusted MH_ratio (MH_ratio_a): %.6f\n', MH_ratio_a);
    fprintf(fid, '\nFinal Host Compartment Values: Short Time\n');
    for j = 1:length(hostData_short_final)
        fprintf(fid, '  x(%d): %.6f\n', j, hostData_short_final(j));
    end
    fprintf(fid, '\nFinal Host Densities: Short Time\n');
    for j = 1:length(total_short_final)
        fprintf(fid, '  H_%d: %.6f\n', j, total_short_final(j));
    end

    fprintf(fid, '\nFinal Host Compartment Values: Long Time\n');
    for j = 1:length(hostData_long_final)
        fprintf(fid, '  x(%d): %.6f\n', j, hostData_long_final(j));
    end
    fprintf(fid, '\nFinal Host Densities: Long Time\n');
    for j = 1:length(total_short_final)
        fprintf(fid, '  H_%d: %.6f\n', j, total_long_final(j));
    end
    fclose(fid);
end

end