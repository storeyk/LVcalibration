% Function to be used with 'main_CA_vary_interactions.m', which simulates
% a CA model described in Cho, Lewis, Storey, Byrne, "Designing 
% experimental conditions to use the Lotka-Volterra model to infer tumor 
% cell line interaction types", 2022.
%
% This function returns the duration of the cell cycle, generated from a
% Gaussian distribution, for a new control or resistant cell in the CA 
% model.
%
% Author: Kathleen Storey <storeyk@lafayette.edu>
% Last revision: 10-18-2022

%%
function cycleDur = initializeCycleDur(state,Cm,Csd,Rm,Rsd)
% Cm,Csd = mean,sd for control cell cycle duration
% Rm,Rsd = mean,sd for resistant cell cycle duration
meanCycleCtrl = Cm;
stdCycleCtrl  = Csd;

meanCycleRes  = Rm;
stdCycleRes   = Rsd;

if state=='R'
    meanCycleDur = meanCycleRes;
    stdCycleDur = stdCycleRes;
elseif state=='C'
    meanCycleDur = meanCycleCtrl;
    stdCycleDur = stdCycleCtrl;
end

cycleDur = -1;
while cycleDur<0
    cycleDur = floor(normrnd(meanCycleDur,stdCycleDur));
end
end