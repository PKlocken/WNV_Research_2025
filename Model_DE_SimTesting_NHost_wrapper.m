clc; clear;

N = 4; % Number of host types
num_years = 5; % Number of years to simulate

urban = [1; 2]; % Urban parameters (1 = less urban, 2 = more urban)
transmit = [0; 1; 34; 4]; % Transmission restriction codes
%origin = ['e'; 'm'; '1'; '2'; '3'; '4']; % Infection origin identifiers
origin = ['m']; % Infection origin identifiers

% Create a timestamped directory to save all figures
timestamp = string(datetime('now', 'Format', 'yyyy_MM_dd_HHmmss'));
baseDir = fullfile('Plots', "Run_" + timestamp);
mkdir(baseDir);

% Loop through all combinations
for j = 1:length(urban)
    for k = 1:length(transmit)
        for i = 1:length(origin)

            orgn = origin(i); % Get current infection origin (char)

            % Create subfolder for this combination
            comboSubDir = fullfile(baseDir, ...
                sprintf('urban_%d__transmission_%d__origin_%s', urban(j), transmit(k), orgn));
            if ~exist(comboSubDir, 'dir')
                mkdir(comboSubDir);
            end

            % Display progress
            fprintf('Running: Urban=%d | Transmission=%d | Origin=%s\n', urban(j), transmit(k), orgn);

            % Call non-seasonal version
            model_DE_SimTesting_NHost_Function(N, urban(j), transmit(k), orgn, num_years, comboSubDir);

            % OPTIONAL: Call seasonal version
            Model_DE_SimTesting_NHost_Function_seasons(N, urban(j), transmit(k), orgn, num_years, comboSubDir);
        end
    end
end

fprintf('\nAll model runs completed.\nFigures saved in directory: %s\n', baseDir);
