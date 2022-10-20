
%This function runs the individual calibration algorithm described in Cho,
%Lewis, Storey, Byrne, "Designing experimental conditions to use the
% Lotka-Volterra model to infer tumor cell line interaction types", 2022.
%
%   Inputs:
%       -type: uses Lotka-Volterra synthetic data 
%               - enter 'Competitive','Mutual', or 'RantagC'
%
%   Note: This script requires use of Marko Laine's MCMC toolbox code
%   (found at https://mjlaine.github.io/mcmcstat/#orgcdeadeb)
%
% Author: Allison Lewis <lewisall@lafayette.edu>
% Last revision: 10-18-2022

function LVcalibration_Individual_LVData( type )


%Inputs:
% type: Competitive, Mutual, RantagC

dataLocation = ['Data/LVData/' type ];
figLocation = ['Figures/Individual_LVData/' type '/'];
mkdir(figLocation)

chainLength = 10000;

global myData weeklyIdx ctrlprop

errorMatrix = zeros(9,1); %keep track of sum of SSEs across all nine data sets
paramsList = zeros(9,7); %keep track of params for each mixture combo
signMatrix = zeros(9,1);
count = 0;

for p1 = 1:9
    
    count = count+1; %tracker for params list
    fprintf('Currently on scenario %d of %d\n', count,9)
    
    
    ctrlprop = 0.1*[p1];
    paramsList(count,1) = ctrlprop; %first entry in param list is the ctrl prop mixture used
    
    
    conds = 0.02.*[ctrlprop 1-ctrlprop]; %get resistant proportions to match ctrl props
    
    
    %Load synthetic LV data
    
    filePattern = fullfile(dataLocation,'*.mat');
    myFiles = dir(filePattern);
    
    
    myFiles2(1) = myFiles(p1);
    
    
    numMixtures = length(myFiles2);
    
    %Loading first patient data
    baseFileName = myFiles2(1).name;
    fullFileName = fullfile(myFiles2(1).folder,baseFileName);
    load(fullFileName)
    
    
    weeklyIdx = [8 15 22 29 36 43 50 57];
    timedata = data.xdata(weeklyIdx)';
    tumordata = data.ydata(:,weeklyIdx)';
    clear data;
    
    
    myData.xdata = timedata;
    myData.ydata(:,1) = tumordata(:,1);
    myData.ydata(:,2) = tumordata(:,2);
    
    %Get initial frequentist estimate
    params = [.5 .5 .5 .5 0 0]; %initial parameter guess
    
    lb = [0 0 0 0 -3 -3];
    ub = [1 1 1 1 3 3];
    
    
    try
        [freqparams, ss0] = fmincon(@sse_lv, params, [], [], [], [], lb, ub)
        
        
        %Run DRAM to get bayesian posteriors
        mse = ss0/(length(myData.ydata(:,1))*2*numMixtures-length(params));
        
        % Create parameter structure
        params1 = {
            {'r_S',freqparams(1),lb(1),ub(1)}
            {'r_R',freqparams(2),lb(2),ub(2)}
            {'K_S', freqparams(3),lb(3),ub(3)}
            {'K_R', freqparams(4),lb(4),ub(4)}
            {'\gammaS',freqparams(5),lb(5),ub(5)}
            {'\gammaR',freqparams(6),lb(6),ub(6)}
            };
        
        data.xdata = myData.xdata;
        data.ydata = myData.ydata;
        
        
        model.ssfun = @sse_lv;
        model.sigma2 = mse;
        options.updatesigma = 1;
        options.method = 'dram';
        no_smps = chainLength;
        
        options.nsimu = no_smps;
        [results,chain,s2chain] = mcmcrun(model,data,params1,options);
        
        options.nsimu = no_smps;
        [results,chain,s2chain,ss2chain] = mcmcrun(model,data,params1,options,results);
        
        
        % Find the optimal parameters
        ind = find(ss2chain == min(ss2chain));  ind = ind(1);
        bayesparams = chain(ind,:) %These are your fitted values
        
        
        resultsFileName = [figLocation 'ChainInfo_MixPt' num2str(p1) '.mat'];
        save(resultsFileName,'results','chain','s2chain','ss2chain','bayesparams')
        
        
        
        
        paramsList(count,2:7) = bayesparams; %add estimated params to list
        
        
        %Did we get the right signs for interaction parameters?
        %Comp: both positive
        %Mut: both negative
        %C antag R: 5 (gammaS) positive, 6 (gammaR) negative
        %R antag C: 5 (gammaS) negative, 6 (gammaR) positive
        switch lower(type)
            case {'competitive'}
                if(bayesparams(5)>0 && bayesparams(6)>0)
                    signMatrix(p1) = 1;
                else
                    signMatrix(p1) = 0;
                end
            case {'mutual'}
                if(bayesparams(5)<0 && bayesparams(6)<0)
                    signMatrix(p1) = 1;
                else
                    signMatrix(p1) = 0;
                end
            case {'cantagr'}
                if(bayesparams(5)>0 && bayesparams(6)<0)
                    signMatrix(p1) = 1;
                else
                    signMatrix(p1) = 0;
                end
            case {'rantagc'}
                if(bayesparams(5)<0 && bayesparams(6)>0)
                    signMatrix(p1) = 1;
                else
                    signMatrix(p1) = 0;
                end
        end
        
        %% Plot graph to compare fit
        figure(1)
        
        % Evaluate model at optimal params
        v0 = conds;
        [time,volume] = ode23(@(t,v)tumorTwoComp(t,v, bayesparams), 0:1:70, v0);
        
        plot(data.xdata,data.ydata(:,1),'ob')
        hold on
        plot(data.xdata,data.ydata(:,2),'or')
        plot(time,volume(:,1),'-b')
        plot(time,volume(:,2),'-r')
        xlabel('Time (days)','FontSize',14);
        ylabel('Volume','FontSize',14);
        set(gca, 'FontSize',14)
        legend('Control','Resistant')
        figName = [figLocation 'Mix_' num2str(p1) '_ModelFits.fig'];
        saveas(gcf, figName)
                
                
        
        figure(2)
        mcmcplot(chain,[],{'r1','r2','K1','K2','\gamma_S','\gamma_R'},'chainpanel')
        figName = [figLocation 'Mix_' num2str(p1) '_Chains.fig'];
        saveas(gcf, figName)
        
        figure(3)
        mcmcplot(chain,[],{'r1','r2','K1','K2','\gamma_S','\gamma_R'},'denspanel')
        figName = [figLocation 'Mix_' num2str(p1) '_Densities.fig'];
        saveas(gcf, figName)
        
        figure(4)
        mcmcplot(chain,[],{'r1','r2','K1','K2','\gamma_S','\gamma_R'},'pairs')
        figName = [figLocation 'Mix_' num2str(p1) '_PairwisePlots.fig'];
        saveas(gcf, figName)
        
        close all
        
        
        %% Check fits for other ratios using these parameters
        conds = .02*[.1 .9; .2 .8; .3 .7; .4 .6; .5 .5; .6 .4; .7 .3; .8 .2; .9 .1];
        
        error = zeros(1,9);
        
        for i = 1:9
            
            %Loading patient data
            baseFileName = myFiles(i).name;
            fullFileName = fullfile(myFiles(i).folder,baseFileName);
            load(fullFileName)
            
            weeklyIdx = [8 15 22 29 36 43 50 57];
            v0 = conds(i,:);
            [time,volume] = ode23(@(t,v)tumorTwoComp(t,v, bayesparams), 0:1:70, v0);
            
            SSE = sum(sum((data.ydata(:,weeklyIdx)'-volume(weeklyIdx,:)).^2));
            error(i) = SSE;
        end
        
        errorMatrix(p1) = sum(error);
        
        
    catch
        errorMatrix(p1) = NaN;
        paramsList(count,2:7) = [NaN NaN NaN NaN NaN NaN];
        signMatrix(p1) = NaN;
    end
    
    
    
    
    
    fileName = [figLocation 'Param&ErrorInfo.mat'];
    save(fileName,'paramsList','errorMatrix','signMatrix')
    
    
end
end


function SSE = sse_lv(params, data)
global myData weeklyIdx ctrlprop

SSE = 0;
conds = 0.02.*[ctrlprop 1-ctrlprop];

v0 = conds;
[time,volume] = ode23(@(t,v)tumorTwoComp(t,v, params), 0:1:70, v0);

if length(volume(:,1))<71
    SSE = 1000; %heavy penalty to get it away from parameter samples that are forcing solution to blow up
else
    SSE = SSE+sum(sum((myData.ydata-volume(weeklyIdx,:)).^2));
end
end
