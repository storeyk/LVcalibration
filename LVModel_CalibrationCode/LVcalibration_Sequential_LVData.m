
%This function runs the sequential algorithm described in Cho,
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

function LVcalibration_Sequential_LVData( type )


global setNum myDataCtrl myDataRes myDataMix fixParams ctrlprop weeklyIdx
no_smps = 10000;

dataLocation = ['Data/LVData/' type ];
figLocation = ['Figures/Sequential_LVData/' type '/'];
mkdir(figLocation)

weeklyIdx = [8 15 22 29 36 43 50 57];


% Control only
setNum = 1;
load('Data/LVData/PureData/CtrlProp_1.mat')


myDataCtrl.xdata = data.xdata(weeklyIdx);
myDataCtrl.ydata = data.ydata(weeklyIdx,:);

%%

params = [.5 .5]; %initial parameter guess

lb = [0 0];
ub = [1 1];
[freqparams, ss0] = fmincon(@sse_lv, params, [], [], [], [], lb, ub)


%Run DRAM to get bayesian posteriors
mse = ss0/(numel(myDataCtrl.ydata)-length(params));

% Create parameter structure
params1 = {
    {'r1',freqparams(1),lb(1),ub(1)}
    {'K1', freqparams(2),lb(2),ub(2)}
    };

data.xdata = myDataCtrl.xdata;
data.ydata = myDataCtrl.ydata;


model.ssfun = @sse_lv;
model.sigma2 = mse;
options.updatesigma = 1;
options.method = 'dram';


options.nsimu = no_smps;
[results,chain,s2chain] = mcmcrun(model,data,params1,options);

options.nsimu = no_smps;
[results,chain,s2chain,ss2chain] = mcmcrun(model,data,params1,options,results);


% Find the optimal parameters
ind = find(ss2chain == min(ss2chain));  ind = ind(1);
bayesparams = chain(ind,:) %These are your fitted values


figure(1)
v0 = [.02 0];
paramVec = [bayesparams(1) 1 bayesparams(2) 1 1 1];
[time,volume] = ode23(@(t,v)tumorTwoComp(t,v, paramVec), 0:1:70, v0);

subplot(2,1,setNum)
plot(data.xdata,data.ydata(:,1),'ob')
hold on
plot(data.xdata,data.ydata(:,2),'or')
plot(time,volume(:,1),'-b')
plot(time,volume(:,2),'-r')
xlabel('Time (days)','FontSize',18);
ylabel('Volume','FontSize',18);
legend('Control','Resistant')


figure(2)
mcmcplot(chain,[],[],'chainpanel')
figName = [figLocation 'PureSensitive_Chains.fig'];
saveas(gcf, figName)

figure(3)
mcmcplot(chain,[],[],'denspanel')
figName = [figLocation 'PureSensitive_Densities.fig'];
saveas(gcf, figName)

figure(4)
mcmcplot(chain,[],[],'pairs')
figName = [figLocation 'PureSensitive_PairwisePlots.fig'];
saveas(gcf, figName)


fixParams = [bayesparams(1) 1 bayesparams(2) 1];


clear data myDataCtrl bayesparams


% NOW CALIBRATE RESISTANT ONLY

setNum = 2;

load('Data/LVData/PureData/CtrlProp_0.mat')

myDataRes.xdata = data.xdata(weeklyIdx);
myDataRes.ydata = data.ydata(weeklyIdx,:);

params = [.5 .5]; %initial parameter guess

lb = [0 0];
ub = [1 1];
[freqparams, ss0] = fmincon(@sse_lv, params, [], [], [], [], lb, ub)


%Run DRAM to get bayesian posteriors
mse = ss0/(numel(myDataRes.ydata)-length(params));

% Create parameter structure
params1 = {
    {'r2',freqparams(1),lb(1),ub(1)}
    {'K2', freqparams(2),lb(2),ub(2)}
    };

data.xdata = myDataRes.xdata;
data.ydata = myDataRes.ydata;


model.ssfun = @sse_lv;
model.sigma2 = mse;
options.updatesigma = 1;
options.method = 'dram';


options.nsimu = no_smps;
[results,chain,s2chain] = mcmcrun(model,data,params1,options);

options.nsimu = no_smps;
[results,chain,s2chain,ss2chain] = mcmcrun(model,data,params1,options,results);


% Find the optimal parameters
ind = find(ss2chain == min(ss2chain));  ind = ind(1);
bayesparams = chain(ind,:) %These are your fitted values


figure(1)
v0 = [0 .02];
paramVec = [1 bayesparams(1) 1 bayesparams(2) 1 1];
[time,volume] = ode23(@(t,v)tumorTwoComp(t,v, paramVec), 0:1:70, v0);

subplot(2,1,setNum)
plot(data.xdata,data.ydata(:,1),'ob')
hold on
plot(data.xdata,data.ydata(:,2),'or')
plot(time,volume(:,1),'-b')
plot(time,volume(:,2),'-r')
xlabel('Time (days)','FontSize',18);
ylabel('Volume','FontSize',18);
legend('Control','Resistant')
figName = [figLocation 'PureExperiments_ModelFits.fig'];
saveas(gcf, figName)

figure(5)
mcmcplot(chain,[],[],'chainpanel')
figName = [figLocation 'PureResistant_Chains.fig'];
saveas(gcf, figName)

figure(6)
mcmcplot(chain,[],[],'denspanel')
figName = [figLocation 'PureResistant_Densities.fig'];
saveas(gcf, figName)


figure(7)
mcmcplot(chain,[],[],'pairs')
figName = [figLocation 'PureResistant_PairwisePlots.fig'];
saveas(gcf, figName)


fixParams(2) = bayesparams(1);
fixParams(4) = bayesparams(2);

clear data myDataRes bayesparams

close all


% NOW CALIBRATE MIXTURE PARAMETERS
setNum = 3;
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
    
    
    timedata = data.xdata(weeklyIdx)';
    tumordata = data.ydata(:,weeklyIdx)';
    clear data;
    
    
    myDataMix.xdata = timedata;
    myDataMix.ydata(:,1) = tumordata(:,1);
    myDataMix.ydata(:,2) = tumordata(:,2);
    
    %Get initial frequentist estimate
    params = [0 0]; %initial parameter guess
    
    lb = [-3 -3];
    ub = [3 3];
    
    
    try
        [freqparams, ss0] = fmincon(@sse_lv, params, [], [], [], [], lb, ub)
        
        
        %Run DRAM to get bayesian posteriors
        mse = ss0/(length(myDataMix.ydata(:,1))*2*numMixtures-length(params));
        
        % Create parameter structure
        params1 = {
            {'\gammaS',freqparams(1),lb(1),ub(1)}
            {'\gammaR',freqparams(2),lb(2),ub(2)}
            };
        
        data.xdata = myDataMix.xdata;
        data.ydata = myDataMix.ydata;
        
        
        model.ssfun = @sse_lv;
        model.sigma2 = mse;
        options.updatesigma = 1;
        options.method = 'dram';
        no_smps = no_smps;
        
        options.nsimu = no_smps;
        [results,chain,s2chain] = mcmcrun(model,data,params1,options);
        
        options.nsimu = no_smps;
        [results,chain,s2chain,ss2chain] = mcmcrun(model,data,params1,options,results);
        
        
        % Find the optimal parameters
        ind = find(ss2chain == min(ss2chain));  ind = ind(1);
        bayesparams = chain(ind,:) %These are your fitted values
        
        
        resultsFileName = [figLocation 'ChainInfo_MixPt' num2str(p1) '.mat'];
        save(resultsFileName,'results','chain','s2chain','ss2chain','bayesparams')
        
        
        
        currParams = [fixParams bayesparams];
        paramsList(count,2:7) = currParams; %add estimated params to list
        
        
        %Did we get the right signs for interaction parameters?
        %Comp: both positive
        %Mut: both negative
        %C antag R: 1 (gammaS) positive, 2 (gammaR) negative
        %R antag C: 1 (gammaS) negative, 2 (gammaR) positive
        switch lower(type)
            case {'competitive'}
                if(bayesparams(1)>0 && bayesparams(2)>0)
                    signMatrix(p1) = 1;
                else
                    signMatrix(p1) = 0;
                end
            case {'mutual'}
                if(bayesparams(1)<0 && bayesparams(2)<0)
                    signMatrix(p1) = 1;
                else
                    signMatrix(p1) = 0;
                end
            case {'cantagr'}
                if(bayesparams(1)>0 && bayesparams(2)<0)
                    signMatrix(p1) = 1;
                else
                    signMatrix(p1) = 0;
                end
            case {'rantagc'}
                if(bayesparams(1)<0 && bayesparams(2)>0)
                    signMatrix(p1) = 1;
                else
                    signMatrix(p1) = 0;
                end
        end
        
        %% Plot graph to compare fit
        figure(8)
        
        % Evaluate model at optimal params
        v0 = conds;
        [time,volume] = ode23(@(t,v)tumorTwoComp(t,v, currParams), 0:1:70, v0);
        
        plot(myDataMix.xdata,myDataMix.ydata(:,1),'ob')
        hold on
        plot(myDataMix.xdata,myDataMix.ydata(:,2),'or')
        plot(time,volume(:,1),'-b')
        plot(time,volume(:,2),'-r')
        xlabel('Time (days)','FontSize',14);
        ylabel('Volume','FontSize',14);
        set(gca, 'FontSize',14)
        legend('Control','Resistant')
        hold off
        figName = [figLocation 'Mix_' num2str(p1) '_ModelFits.fig'];
        saveas(gcf, figName)

        
        figure(9)
        mcmcplot(chain,[],{'\gamma_S','\gamma_R'},'chainpanel')
        figName = [figLocation 'Mix_' num2str(p1) '_Chains.fig'];
        saveas(gcf, figName)
        
        figure(10)
        mcmcplot(chain,[],{'\gamma_S','\gamma_R'},'denspanel')
        figName = [figLocation 'Mix_' num2str(p1) '_Densities.fig'];
        saveas(gcf, figName)
        
        figure(11)
        mcmcplot(chain,[],{'\gamma_S','\gamma_R'},'pairs')
        figName = [figLocation 'Mix_' num2str(p1) '_PairwisePlots.fig'];
        saveas(gcf, figName)
        
        
        
        %% Check fits for other ratios using these parameters
        conds = .02*[.1 .9; .2 .8; .3 .7; .4 .6; .5 .5; .6 .4; .7 .3; .8 .2; .9 .1];
        
        error = zeros(1,9);
        
        for i = 1:9
            
            %Loading patient data
            baseFileName = myFiles(i).name;
            fullFileName = fullfile(myFiles(i).folder,baseFileName);
            load(fullFileName)
            
            v0 = conds(i,:);
            [time,volume] = ode23(@(t,v)tumorTwoComp(t,v, currParams), 0:1:70, v0);
            
            SSE = sum(sum((data.ydata(:,weeklyIdx)'-volume(weeklyIdx,:)).^2));
            error(i) = SSE;
        end
        
        errorMatrix(p1) = sum(error);
        
        
    catch
        errorMatrix(p1) = NaN;
        paramsList(count,2:7) = [NaN NaN NaN NaN NaN NaN];
        signMatrix(p1) = NaN;
    end
    
    
    
    
    
    clear data myDataMix bayesparams
end


fileName = [figLocation 'Param&ErrorInfo.mat'];
save(fileName,'paramsList','errorMatrix','signMatrix')


end


function SSE = sse_lv(params, data)
global myDataCtrl myDataRes myDataMix setNum fixParams ctrlprop weeklyIdx

SSE = 0;

if setNum == 1
    v0 = [.02 0];
    paramVec = [params(1) 1 params(2) 1 1 1];
    [time,volume] = ode23(@(t,v)tumorTwoComp(t,v, paramVec), 0:1:70, v0);
    if length(volume(:,1))<71
        SSE = 1000; %heavy penalty to get it away from parameter samples that are forcing solution to blow up
    else
        SSE = sum((myDataCtrl.ydata(:,1)-volume(weeklyIdx,1)).^2);
    end
end
if setNum == 2
    v0 = [0 .02];
    paramVec = [1 params(1) 1 params(2) 1 1];
    [time,volume] = ode23(@(t,v)tumorTwoComp(t,v, paramVec), 0:1:70, v0);
    if length(volume(:,1))<71
        SSE = 1000; %heavy penalty to get it away from parameter samples that are forcing solution to blow up
    else
        SSE = sum((myDataRes.ydata(:,2)-volume(weeklyIdx,2)).^2);
    end
end
if setNum == 3
    v0 = 0.02.*[ctrlprop 1-ctrlprop];
    paramVec = [fixParams params];
    [time,volume] = ode23(@(t,v)tumorTwoComp(t,v, paramVec), 0:1:70, v0);
    if length(volume(:,1))<71
        SSE = 1000; %heavy penalty to get it away from parameter samples that are forcing solution to blow up
    else
        SSE = sum(sum((myDataMix.ydata-volume(weeklyIdx,:)).^2));
    end
end

end
