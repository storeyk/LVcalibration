% This is the main script to simulate the CA model described in Cho,
% Lewis, Storey, Byrne, "Designing experimental conditions to use the
% Lotka-Volterra model to infer tumor cell line interaction types", 2022.
%
%   Outputs:
%       - Matrix of volumes of each cell type, for each day.
%       - Matrix of cell types at each spatial location, for each day.
%
% Author: Kathleen Storey <storeyk@lafayette.edu>
% Last revision: 10-18-2022

%% Create and initialize the states
clear variables; close all; clc;
CELLSIZE       = 0.0018; %[cm]
totTime        = 24*70;
totRep         = 1;
n              = 200;    %grid size


%%Cell cycle times
Cm = 18.3;      %mean Control cell cycle time (PC3 baseline: 18.3)
Csd = 1.4;      %sd Control cell cycle time (PC3 baseline: 1.4)
Rm = 18.3;      %mean Resistant cell cycle time (PC3 baseline: 16.9)
Rsd = 1.4;      %sd Resistant cell cycle time (PC3 baseline: 0.9)

cQ        = 0.8;                   %oxygen threshold for control proliferating 
cq        = 0.8;                   %oxygen threshold for resistant proliferating
cN        = 0.75;                  %oxygen threshold for control quiescent 
cn        = 0.75;                  %oxygen threshold for resistant quiescent
kappaC    = 1e-8;                  %oxygen consumption rate control 
kappaR    = 1e-8;                  %oxygen consumption rate resistant


% list of initial proportions of Type-S cells to loop through:
prop_list = 0.1:0.1:0.9;

for propInd=1:length(prop_list)
ctrl_res_split = prop_list(propInd); % Initial proportion of Type-S cells

%%% Set rad on(1) or off(0):
rad=0;

%%% Set interaction type:
%%0=neutral, 1=competition, 2=mutualism, 3=antagon (C helps R, R hurts C - label: RantagC),
%%4=antagon (R helps C, C hurts R - label:CantagR)
interact=1;

% Set a range of parameter values to loop through 
% Interaction intensity (Iintens) vals:
IIvals=[2,3,4];

% Initialize matrix to store volumes
max_time=totTime/24+1;
vol=zeros((max_time-1),10,totRep);

dir_name=['Interact_data_twoTypes/CtrlProp_' num2str(ctrl_res_split)];

if rad
    if ctrl_res_split == 0 || ctrl_res_split == 1
        
        dir_name = [dir_name '/rad'];
        
    else
         
        if interact==1
            dir_name = [dir_name '/comp'];
        elseif interact==2
            dir_name = [dir_name '/mutual'];
        elseif interact==3
            dir_name = [dir_name '/RantagC'];
        elseif interact==4
            dir_name = [dir_name '/CantagR'];
        else
            dir_name = [dir_name '/neutral'];
        end
    end
else
    if ctrl_res_split == 0 || ctrl_res_split == 1
        
        dir_name = [dir_name '/norad'];
        
    else
    
        if interact==1
            dir_name = [dir_name '/norad_comp'];
        elseif interact==2
            dir_name = [dir_name '/norad_mutual'];
        elseif interact==3
            dir_name = [dir_name '/norad_RantagC'];
        elseif interact==4
            dir_name = [dir_name '/norad_CantagR'];
        else
            dir_name = [dir_name '/norad_neutral'];
        end
    end
end
mkdir(dir_name)


for pind=1:length(IIvals)

% interaction parameters
Thresh1 = 4;%4        % Threshold 1: number of total nbrs over which the div ctr reduces by tau/2 in neutral case
                      % and number of nbrs of opposite type affects reduction of div counter
Iintens = IIvals(pind)      % controls the interaction intensity (level of reduction of div counter, in comparison to neutral case)
                            % neutral case: Iintens=1 

% Related (but not specifically interaction) params
Thresh2 = 7;        % Threshold 2: number of total nbrs over which the div ctr reduces by tau/4 in neutral case
                    % (no change for interaction effect), we need Thresh1<Thresh2
                    
NbrhdSize = 8;      % Number of neighboring sites to check (starting with 8) 


                                     
filename = [dir_name '/T1_' num2str(Thresh1) '_I_' num2str(Iintens)];



    p_NR      = 0.0075;%pNRvals(pind);
    %Assuming we want same necrotic decay rate for both:
    p_nR      = 0.0075;%pNRvals(pind); 
for rep=1:totRep
    %rep
    %% Set and scale parameter values
    tau       = 0.5;                   %time step for CA
    p_PQ      = 0.25;                  %probability of transition from P to Q
    p_QP      = 0.25;                  %probability of transition from Q back to P
    p_QN      = 0.25;                  %probability of transition from Q to N
    
    
    
    %----------------------------------------------------------------------
    % These parameter values can be edited
    if rad
        dose = 2;                      %rad dose per fraction
    else
        dose = 0;
    end
    alphaC    = 0.14;                   %radiosensitivity parameter alpha for control
    alphaR    = 0.14;                   %radiosensitivity parameter alpha for ressistant
    betaC     = 0.0467;                 %radiosensitivity parameter beta for control
    betaR     = 0.0156;                 %radiosensitivity parameter beta for resistant
    
    
    numDoses  = 30;                    %total number of rad doses
    firstDose = 15*24+1;               %time first dose is given
    %-----------------------------------------------------------------------
    
    L         = n*CELLSIZE;            %domain length [cm]
    d         = 1;                     %oxygen diffusion rate
    delx      = CELLSIZE/L;            %RDE space step or 1/n
    delt      = delx^2/(4*d);          %RDE time step
    tol       = 1.0e-6;                %FD method convergence tolerance
    courant   = (d*delt)/delx^2;       %Courant number (must be <= 0.5)
    cInf      = 2.8e-7;                %background oxygen concentration [mol/cm^3]
    bc        = 1;                     %O2 boundary conditions
    oxygen    = ones(n);               %initial oxygen concentration
    D         = 1.8e-5;                %oxygen diffusion constant [cm^2/s]
    kC        = (L^2*kappaC)/(cInf*D); %prolif cells O2 consumption rate control
    kR        = (L^2*kappaR)/(cInf*D); %prolif cells O2 consumption rate resistant
    kQ        = 0.5*kC;                %quiescent cells O2 consumption rate control
    kq        = 0.5*kR;                %quiescent cells O2 consumption rate resistant
    
    p_radC    = 1-exp(-alphaC*dose-betaC*dose^2);   %probability of proliferating control cell death via radiation
    p_radR    = 1-exp(-alphaR*dose-betaR*dose^2);   %probability of proliferating resistant cell death via radiation
    
    
    %% Initialize grid with cells
    r     = 9;
    h     = n/2;  %center of the circle
    state = char(zeros(n,n));
    for x=1:n
        for y=1:n
            if (x-h)^2+(y-h)^2<=r^2
                if rand<ctrl_res_split
                    state(x,y)='C';
                else
                    state(x,y)='R';
                end
            else
                state(x,y)='E';
            end
        end
    end
    % Randomly initialize cell cycle times
    divisionCounter = zeros(n);
    divisionCounter(state=='C') = abs(floor(18.3.*rand(size(divisionCounter(state=='C')))));
    divisionCounter(state=='R') = abs(floor(16.9.*rand(size(divisionCounter(state=='R')))));
    
    %% Store states at different time points
    t=1;
    allStates = char(zeros(n,n,totTime/24+1));
    allStates(:,:,t) = state;
    t=t+1;
    
    %% Initialize rad doses
    nextDose = firstDose;       %time of next dose
    curDose = 0;                %number of doses already given 
    
    %% Update grid
    for iTime=0.5:tau:totTime
        
        tmp = zeros(1,n);
        tmp(state(n/2,:)~='E') = 1;
        if sum(tmp)>=180
            break
        end
        
        %% Update oxygen concentration
        % Update consumption rates based on cell state
        consumptionRate             = zeros(n);
        consumptionRate(state=='C') = kC;
        consumptionRate(state=='R') = kR;
        consumptionRate(state=='Q') = kQ;
        consumptionRate(state=='q') = kq;
        
        err = 1;
        totalT=0;
        tmp = zeros(n);
        while err > tol
            totalT=totalT+delt;
            tmp(1,:) = bc;
            tmp(:,1) = bc;
            tmp(n,:) = bc;
            tmp(:,n) = bc;
            for i = 2:n-1
                for j = 2:n-1
                    if state(i,j)=='E' %BCs
                        tmp(i,j)=bc;
                    else
                        tmp(i,j) = oxygen(i,j)+courant*(oxygen(i-1,j)+oxygen(i+1,j)+oxygen(i,j-1)+oxygen(i,j+1)-4*oxygen(i,j))-delt*consumptionRate(i,j);
                    end
                end
            end
            err = max(max(abs(tmp-oxygen)));
            oxygen = tmp;
        end
        
        %% Update cell states based on oxygen concentration
        % If O2 level below threshold change state to Q
        tmp = ones(n);
        tmp(state=='C' & oxygen<cQ) = rand(size(tmp(state=='C' & oxygen<cQ)));
        state(tmp<p_PQ) = 'Q';
        
        tmp = ones(n);
        tmp(state=='R' & oxygen<cq) = rand(size(tmp(state=='R' & oxygen<cq)));
        state(tmp<p_PQ) = 'q';
        
        % If O2 level above threshold change state to P
        tmp = ones(n);
        tmp(state=='Q' & oxygen>cQ) = rand(size(tmp(state=='Q' & oxygen>cQ)));
        state(tmp<p_QP) = 'C';
        
        tmp = ones(n);
        tmp(state=='q' & oxygen>cq) = rand(size(tmp(state=='q' & oxygen>cq)));
        state(tmp<p_QP) = 'R';
        
        % If O2 level below threshold change state to N
        tmp = ones(n);
        tmp(state=='Q' & oxygen<cN) = rand(size(tmp(state=='Q' & oxygen<cN)));
        state(tmp<p_QN) = 'N';
        divisionCounter(tmp<p_QN) = 0;
        
        tmp = ones(n);
        tmp(state=='q' & oxygen<cn) = rand(size(tmp(state=='q' & oxygen<cn)));
        state(tmp<p_QN) = 'n';
        divisionCounter(tmp<p_QN) = 0;
        
        %% Radiation dose 
        if (curDose<numDoses) && (iTime > nextDose)
            tmp = ones(n);
            tmp(state=='C') = rand(size(tmp(state=='C')));
            state(tmp<p_radC) = 'N';
            
            tmp = ones(n);
            tmp(state=='R') = rand(size(tmp(state=='R')));
            state(tmp<p_radR) = 'n';
            
            %quiescent cells: (2/3)*prob of rad death 
            tmp = ones(n);
            tmp(state=='Q') = rand(size(tmp(state=='Q')));
            state(tmp<(2/3)*p_radC) = 'N';
            
            tmp = ones(n);
            tmp(state=='q') = rand(size(tmp(state=='q')));
            state(tmp<(2/3)*p_radR) = 'n';
            
            
            curDose = curDose+1;
            if mod(curDose,5)==0
               nextDose = nextDose+72;  
            else
               nextDose = nextDose+24; 
            end
        end
                
        
        %% Store results
        if mod(iTime,24)==0
            allStates(:,:,t) = state;
            t=t+1
        end
        
        %% Update cell cycle time
        % Find out what phase my neighbors are in
        [totalFree, totalOpp] = checkNeighbors(state,n);
        totalNbrs = NbrhdSize - totalFree;
        
        % Reduce division counter
        divisionCounter((state=='C'|state=='R') & totalNbrs<=Thresh1) = divisionCounter((state=='C'|state=='R') & totalNbrs<=Thresh1)-tau/1;
        %cell consumes oxygen proportional to division reduction
        oxygen((state=='C'|state=='R') & totalNbrs<=Thresh1) = oxygen((state=='C'|state=='R') & totalNbrs<=Thresh1) - delt*kC;
        
        % Between thresh1 and thresh2, but no interaction effect (number of nbrs of opp type is <=4)
        divisionCounter((state=='C'|state=='R') & totalNbrs>Thresh1 & totalNbrs<=Thresh2 & totalOpp<=Thresh1) = divisionCounter((state=='C'|state=='R') & totalNbrs>Thresh1 & totalNbrs<=Thresh2 & totalOpp<=Thresh1)-tau/2;
        oxygen((state=='C'|state=='R') & totalNbrs<=Thresh1 & totalNbrs<=Thresh2 & totalOpp<=Thresh1) = oxygen((state=='C'|state=='R') & totalNbrs<=Thresh1 & totalNbrs<=Thresh2 & totalOpp<=Thresh1) - delt/2*kC;
        % Above thresh2, no interaction effect
        divisionCounter((state=='C'|state=='R') & totalNbrs>Thresh2 & totalOpp<=Thresh1) = divisionCounter((state=='C'|state=='R') & totalNbrs>Thresh2 & totalOpp<=Thresh1)-tau/4;
        oxygen((state=='C'|state=='R') & totalNbrs>Thresh2 & totalOpp<=Thresh1) = oxygen((state=='C'|state=='R') & totalNbrs>Thresh2 & totalOpp<=Thresh1) - delt/4*kC;
        
         
        if interact==1
            %divide cell counter reduction by Iintens
            divisionCounter((state=='C'|state=='R') & totalOpp>Thresh1 & totalNbrs<=Thresh2) = divisionCounter((state=='C'|state=='R') & totalOpp>Thresh1 & totalNbrs<=Thresh2)-tau/(2*Iintens);
            divisionCounter((state=='C'|state=='R') & totalOpp>Thresh1 & totalNbrs>Thresh2) = divisionCounter((state=='C'|state=='R') & totalOpp>Thresh1 & totalNbrs>Thresh2)-tau/(4*Iintens);
            oxygen((state=='C'|state=='R') & totalOpp>Thresh1 & totalNbrs<=Thresh2) = oxygen((state=='C'|state=='R') & totalOpp>Thresh1 & totalNbrs<=Thresh2)-delt/(2*Iintens)*kC;
            oxygen((state=='C'|state=='R') & totalOpp>Thresh1 & totalNbrs>Thresh2) = oxygen((state=='C'|state=='R') & totalOpp>Thresh1 & totalNbrs>Thresh2)-delt/(4*Iintens)*kC;
            
         elseif interact==2
            %multiply cell counter by Iintens
            divisionCounter((state=='C'|state=='R') & totalOpp>Thresh1 & totalNbrs<=Thresh2) = divisionCounter((state=='C'|state=='R') & totalOpp>Thresh1 & totalNbrs<=Thresh2)-tau/2*Iintens;
            divisionCounter((state=='C'|state=='R') & totalOpp>Thresh1 & totalNbrs>Thresh2) = divisionCounter((state=='C'|state=='R') & totalOpp>Thresh1 & totalNbrs>Thresh2)-tau/4*Iintens;
            oxygen((state=='C'|state=='R') & totalOpp>Thresh1 & totalNbrs<=Thresh2) = oxygen((state=='C'|state=='R') & totalOpp>Thresh1 & totalNbrs<=Thresh2)-delt/2*Iintens*kC;
            oxygen((state=='C'|state=='R') & totalOpp>Thresh1 & totalNbrs>Thresh2) = oxygen((state=='C'|state=='R') & totalOpp>Thresh1 & totalNbrs>Thresh2)-delt/4*Iintens*kC;
            
            
        elseif interact==3
            %C helps R, R hurts C
            divisionCounter(state=='C' & totalOpp>Thresh1 & totalNbrs<=Thresh2) = divisionCounter(state=='C' & totalOpp>Thresh1 & totalNbrs<=Thresh2)-tau/(2*Iintens);
            divisionCounter(state=='C' & totalOpp>Thresh1 & totalNbrs>Thresh2) = divisionCounter(state=='C' & totalOpp>Thresh1 & totalNbrs>Thresh2)-tau/(4*Iintens);
            oxygen(state=='C' & totalOpp>Thresh1 & totalNbrs<=Thresh2) = oxygen(state=='C' & totalOpp>Thresh1 & totalNbrs<=Thresh2)-delt/(2*Iintens)*kC;
            oxygen(state=='C' & totalOpp>Thresh1 & totalNbrs>Thresh2) = oxygen(state=='C' & totalOpp>Thresh1 & totalNbrs>Thresh2)-delt/(4*Iintens)*kC;
            
            divisionCounter(state=='R' & totalOpp>Thresh1 & totalNbrs<=Thresh2) = divisionCounter(state=='R' & totalOpp>Thresh1 & totalNbrs<=Thresh2)-tau/2*Iintens;
            divisionCounter(state=='R' & totalOpp>Thresh1 & totalNbrs>Thresh2) = divisionCounter(state=='R' & totalOpp>Thresh1 & totalNbrs>Thresh2)-tau/4*Iintens;
            oxygen(state=='R' & totalOpp>Thresh1 & totalNbrs<=Thresh2) = oxygen(state=='R' & totalOpp>Thresh1 & totalNbrs<=Thresh2)-delt/2*Iintens*kC;
            oxygen(state=='R' & totalOpp>Thresh1 & totalNbrs>Thresh2) = oxygen(state=='R' & totalOpp>Thresh1 & totalNbrs>Thresh2)-delt/4*Iintens*kC;
            
            
        elseif interact==4
            %C hurts R, R helps C
            divisionCounter(state=='R' & totalOpp>Thresh1 & totalNbrs<=Thresh2) = divisionCounter(state=='R' & totalOpp>Thresh1 & totalNbrs<=Thresh2)-tau/(2*Iintens);
            divisionCounter(state=='R' & totalOpp>Thresh1 & totalNbrs>Thresh2) = divisionCounter(state=='R' & totalOpp>Thresh1 & totalNbrs>Thresh2)-tau/(4*Iintens);
            oxygen(state=='R' & totalOpp>Thresh1 & totalNbrs<=Thresh2) = oxygen(state=='R' & totalOpp>Thresh1 & totalNbrs<=Thresh2)-delt/(2*Iintens)*kC;
            oxygen(state=='R' & totalOpp>Thresh1 & totalNbrs>Thresh2) = oxygen(state=='R' & totalOpp>Thresh1 & totalNbrs>Thresh2)-delt/(4*Iintens)*kC;
            
            divisionCounter(state=='C' & totalOpp>Thresh1 & totalNbrs<=Thresh2) = divisionCounter(state=='C' & totalOpp>Thresh1 & totalNbrs<=Thresh2)-tau/2*Iintens;
            divisionCounter(state=='C' & totalOpp>Thresh1 & totalNbrs>Thresh2) = divisionCounter(state=='C' & totalOpp>Thresh1 & totalNbrs>Thresh2)-tau/4*Iintens;
            oxygen(state=='C' & totalOpp>Thresh1 & totalNbrs<=Thresh2) = oxygen(state=='C' & totalOpp>Thresh1 & totalNbrs<=Thresh2)-delt/2*Iintens*kC;
            oxygen(state=='C' & totalOpp>Thresh1 & totalNbrs>Thresh2) = oxygen(state=='C' & totalOpp>Thresh1 & totalNbrs>Thresh2)-delt/4*Iintens*kC;
            
            
        end
        
        %% Remove necrotic cells
        % Flag cells to be removed with a rand
        tmp = ones(n);
        tmp(state=='N') = rand(size(tmp(state=='N')));
        removeFlag = zeros(n);
        removeFlag(tmp<p_NR) = rand(size(removeFlag(tmp<p_NR)));
        
        % Iterate over cells being removed starting from the largerst rand
        while nnz(removeFlag)~=0
            [row,col] = find(removeFlag==max(max(removeFlag))); %indices of cell with highest rand
            removeFlag(row,col)=0;
            [state,divisionCounter,removeFlag] = removeDead(state,divisionCounter,removeFlag,row,col,n);
        end
        
        % Flag cells to be removed with a rand
        tmp = ones(n);
        tmp(state=='n') = rand(size(tmp(state=='n')));
        removeFlag = zeros(n);
        removeFlag(tmp<p_nR) = rand(size(removeFlag(tmp<p_nR)));
        
        % Iterate over cells being removed starting from the largerst rand
        while nnz(removeFlag)~=0
            [row,col] = find(removeFlag==max(max(removeFlag))); %indices of cell with highest rand
            removeFlag(row,col)=0;
            [state,divisionCounter,removeFlag] = removeDead(state,divisionCounter,removeFlag,row,col,n);
        end
        
        %% Actions performed by proliferating cells
        % If division counter is zero then attempt division
        divisionFlag = zeros(n);
        divisionFlag((state=='C'|state=='R') & divisionCounter<=0) = rand(size(divisionFlag((state=='C'|state=='R') & divisionCounter<=0)));
        
        % Iterate over dividing cells starting from the largerst rand
        while nnz(divisionFlag)~=0
            [row,col] = find(divisionFlag==max(max(divisionFlag))); %indices of cell with highest rand
            divisionFlag(row,col)=0;
            divisionCounter(row,col)=initializeCycleDur(state(row,col),Cm,Csd,Rm,Rsd);
            [state,divisionCounter,divisionFlag] = divide(state,divisionCounter,divisionFlag,row,col,n,Cm,Csd,Rm,Rsd);
        end
            

    end
    
    for time=1:max_time
        vol(time,:,rep)=[time-1 computeVolume(allStates(:,:,time),n)];
    end
    
end

xdata=0:max_time-1;
ydata=vol(:,:,:);

data.xdata=xdata;
data.ydata=ydata;
writefile=[filename '_vols.mat'];
save(writefile,'data')

heterogeneous=allStates;
fileStates=[filename '_spatial.mat'];
save(fileStates,'heterogeneous')


end

end
