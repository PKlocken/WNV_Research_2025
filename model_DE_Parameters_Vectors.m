% Updated 01/12/2025
% For citations regarding vector parameters, please refer to comments in West Nile code written by Dr. Wandi Ding and Dr. Rachel Leander.

function [pV] = model_DE_Parameters_Vectors()

% Parameter values for simulations:
pV = zeros(1,13);

% VECTOR PARAMETERS

pV(1) = 150/(2*8); % egg laying rate of S and E mosquitoes, rs, 
% assumes 150 eggs per raft, half of these are female, and a gonotrophic
% cycle of 8 days
pV(2) = 100/(2*8); % egg laying rate of I mosquitoes, ri

pV(3) = 0.003; % fraction of eggs infected, phi

pV(4) = .56; % fraction of eggs laid by uninfected mosquitoes that hatch, qs
pV(5) = .43; % fraction of eggs laid by infected mosquitoes that hatch, qi

pV(6) = 1/2; % hatch rate, m_e

pV(7) = 1/7; % larval maturation rate (1/larval lifespan), m_L
            % Here we combine the larval and pupal stages. 

pV(8) = 0.16; % Daily death (of larvae) rate of 1 - 3%, mu_L             
             % Larval stage survival varied between 70-90%. 
             % This would mean p(7)/(p(7)+p(8))=.7-.9.
             % For p(8) and p(7) as above this ratio is .82.                  
                  
pV(9) = 1/10.4; % adult death rate (1/adult lifespan), mu_V
               % female mosquitoes have a life expectancy of 3-7 days in the wild

pV(10) = 1/5; % mosquito biting rate, b

%pV(11) = .01; % mosquito larval carrying capacity, c_L per m^2
pV(11) = 100; % mosquito larval carrying capacity, c_L per ha
%pV(11) = 10; % mosquito larval carrying capacity, c_L per ha

pV(12) = 1/10; % disease progression in mosquitoes (1/latency period), kl

%pV(13) = 0.5; % mosquito-to-host transmission, p_mh
pV(13) = 1; % mosquito-to-host transmission, p_mh

end