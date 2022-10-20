% Function to be used with 'main_CA_vary_interactions.m', which simulates
% a CA model described in Cho, Lewis, Storey, Byrne, "Designing 
% experimental conditions to use the Lotka-Volterra model to infer tumor 
% cell line interaction types", 2022.
%
% This function produces a new daughter cell in a neighboring site, if 
% space is available. If there is no space available, we find the smallest
% distance to the spheroid surface and shift cells in that
% distance by one.
%
% Author: Kathleen Storey <storeyk@lafayette.edu>
% Last revision: 10-18-2022

%%
function [state,divisionCounter,divisionFlag] = divide(state,divisionCounter,divisionFlag,row,col,n,Cm,Csd,Rm,Rsd)

currState=state(row,col); currDivisionCounter=initializeCycleDur(state(row,col),Cm,Csd,Rm,Rsd); currDivisionFlag=0;

while (state(row-1,col)~='E' && state(row-1,col-1)~='E' && state(row-1,col+1)~='E' &&...
        state(row,col-1)~='E' && state(row,col+1)~='E' &&...
        state(row+1,col-1)~='E' && state(row+1,col)~='E' && state(row+1,col+1)~='E')
    
    % Look in all 8 directions and find the direction with lowest distance
    % to the spheroid surface.
    % North...
    i=row; j=col; northDistance=0;
    while state(i-1,j)~='E'
        i=i-1; northDistance=northDistance+1;
    end
    % East...
    i=row; j=col; eastDistance=0;
    while state(i,j+1)~='E'
        j=j+1; eastDistance=eastDistance+1;
    end
    % South...
    i=row; j=col; southDistance=0;
    while state(i+1,j)~='E'
        i=i+1; southDistance=southDistance+1;
    end
    % West...
    i=row; j=col; westDistance=0;
    while state(i,j-1)~='E'
        j=j-1; westDistance=westDistance+1;
    end
    % North-east...
    i=row; j=col; northEastDistance=0;
    while state(i-1,j+1)~='E'
        i=i-1; j=j+1; northEastDistance=northEastDistance+1;
    end
    % North-west...
    i=row; j=col; northWestDistance=0;
    while state(i-1,j-1)~='E'
        i=i-1; j=j-1; northWestDistance=northWestDistance+1;
    end
    % South-east
    i=row; j=col; southEastDistance=0;
    while state(i+1,j+1)~='E'
        i=i+1; j=j+1; southEastDistance=southEastDistance+1;
    end
    % South-west...
    i=row; j=col; southWestDistance=0;
    while state(i+1,j-1)~='E'
        i=i+1; j=j-1; southWestDistance=southWestDistance+1;
    end
    
    % Push cell in the biased random direction
    distances=[northDistance eastDistance southDistance westDistance...
        northEastDistance northWestDistance southEastDistance southWestDistance];
    tmp=zeros(1,8);
    sortedDistances=sort(distances);
    
    tmp(distances==sortedDistances(8))=rand(size(tmp(distances==sortedDistances(8))));
    tmp(distances==sortedDistances(7))=rand(size(tmp(distances==sortedDistances(7))));
    tmp(distances==sortedDistances(6))=rand(size(tmp(distances==sortedDistances(6))));
    tmp(distances==sortedDistances(5))=rand(size(tmp(distances==sortedDistances(5))))+0.6;
    tmp(distances==sortedDistances(4))=rand(size(tmp(distances==sortedDistances(4))))+0.7;
    tmp(distances==sortedDistances(3))=rand(size(tmp(distances==sortedDistances(3))))+0.8;
    tmp(distances==sortedDistances(2))=rand(size(tmp(distances==sortedDistances(2))))+0.9;
    tmp(distances==sortedDistances(1))=rand(size(tmp(distances==sortedDistances(1))))+1;
    
    northDistance=tmp(1);
    eastDistance=tmp(2);
    southDistance=tmp(3);
    westDistance=tmp(4);
    northEastDistance=tmp(5);
    northWestDistance=tmp(6);
    southEastDistance=tmp(7);
    southWestDistance=tmp(8);
    
    %Shift cell along the lowest distance to nearest empty cell
    if northDistance==max([northDistance eastDistance southDistance westDistance...
            northEastDistance northWestDistance southEastDistance southWestDistance])
        tmpState                   = state(row-1,col);
        tmpDivisionCounter         = divisionCounter(row-1,col);
        tmpDivisionFlag            = divisionFlag(row-1,col);
        state(row-1,col)           = currState;
        divisionCounter(row-1,col) = currDivisionCounter;
        divisionFlag(row-1,col)    = currDivisionFlag;
        currState                  = tmpState;
        currDivisionCounter        = tmpDivisionCounter;
        currDivisionFlag           = tmpDivisionFlag;
        row                        = row-1;
    elseif eastDistance==max([northDistance eastDistance southDistance westDistance...
            northEastDistance northWestDistance southEastDistance southWestDistance])
        tmpState                   = state(row,col+1);
        tmpDivisionCounter         = divisionCounter(row,col+1);
        tmpDivisionFlag            = divisionFlag(row,col+1);
        state(row,col+1)           = currState;
        divisionCounter(row,col+1) = currDivisionCounter;
        divisionFlag(row,col+1)    = currDivisionFlag;
        currState                  = tmpState;
        currDivisionCounter        = tmpDivisionCounter;
        currDivisionFlag           = tmpDivisionFlag;
        col                        = col+1;
    elseif southDistance==max([northDistance eastDistance southDistance westDistance...
            northEastDistance northWestDistance southEastDistance southWestDistance])
        tmpState                   = state(row+1,col);
        tmpDivisionCounter         = divisionCounter(row+1,col);
        tmpDivisionFlag            = divisionFlag(row+1,col);
        state(row+1,col)           = currState;
        divisionCounter(row+1,col) = currDivisionCounter;
        divisionFlag(row+1,col)    = currDivisionFlag;
        currState                  = tmpState;
        currDivisionCounter        = tmpDivisionCounter;
        currDivisionFlag           = tmpDivisionFlag;
        row                        = row+1;
    elseif westDistance==max([northDistance eastDistance southDistance westDistance...
            northEastDistance northWestDistance southEastDistance southWestDistance])
        tmpState                   = state(row,col-1);
        tmpDivisionCounter         = divisionCounter(row,col-1);
        tmpDivisionFlag            = divisionFlag(row,col-1);
        state(row,col-1)           = currState;
        divisionCounter(row,col-1) = currDivisionCounter;
        divisionFlag(row,col-1)    = currDivisionFlag;
        currState                  = tmpState;
        currDivisionCounter        = tmpDivisionCounter;
        currDivisionFlag           = tmpDivisionFlag;
        col                        = col-1;
    elseif northEastDistance==max([northDistance eastDistance southDistance westDistance...
            northEastDistance northWestDistance southEastDistance southWestDistance])
        tmpState                     = state(row-1,col+1);
        tmpDivisionCounter           = divisionCounter(row-1,col+1);
        tmpDivisionFlag              = divisionFlag(row-1,col+1);
        state(row-1,col+1)           = currState;
        divisionCounter(row-1,col+1) = currDivisionCounter;
        divisionFlag(row-1,col+1)    = currDivisionFlag;
        currState                    = tmpState;
        currDivisionCounter          = tmpDivisionCounter;
        currDivisionFlag             = tmpDivisionFlag;
        row                          = row-1;
        col                          = col+1;
    elseif northWestDistance==max([northDistance eastDistance southDistance westDistance...
            northEastDistance northWestDistance southEastDistance southWestDistance])
        tmpState                     = state(row-1,col-1);
        tmpDivisionCounter           = divisionCounter(row-1,col-1);
        tmpDivisionFlag              = divisionFlag(row-1,col-1);
        state(row-1,col-1)           = currState;
        divisionCounter(row-1,col-1) = currDivisionCounter;
        divisionFlag(row-1,col-1)    = currDivisionFlag;
        currState                    = tmpState;
        currDivisionCounter          = tmpDivisionCounter;
        currDivisionFlag             = tmpDivisionFlag;
        row                          = row-1;
        col                          = col-1;
    elseif southEastDistance==max([northDistance eastDistance southDistance westDistance...
            northEastDistance northWestDistance southEastDistance southWestDistance])
        tmpState                     = state(row+1,col+1);
        tmpDivisionCounter           = divisionCounter(row+1,col+1);
        tmpDivisionFlag              = divisionFlag(row+1,col+1);
        state(row+1,col+1)           = currState;
        divisionCounter(row+1,col+1) = currDivisionCounter;
        divisionFlag(row+1,col+1)    = currDivisionFlag;
        currState                    = tmpState;
        currDivisionCounter          = tmpDivisionCounter;
        currDivisionFlag             = tmpDivisionFlag;
        row                          = row+1;
        col                          = col+1;
    elseif southWestDistance==max([northDistance eastDistance southDistance westDistance...
            northEastDistance northWestDistance southEastDistance southWestDistance])
        tmpState                     = state(row+1,col-1);
        tmpDivisionCounter           = divisionCounter(row+1,col-1);
        tmpDivisionFlag              = divisionFlag(row+1,col-1);
        state(row+1,col-1)           = currState;
        divisionCounter(row+1,col-1) = currDivisionCounter;
        divisionFlag(row+1,col-1)    = currDivisionFlag;
        currState                    = tmpState;
        currDivisionCounter          = tmpDivisionCounter;
        currDivisionFlag             = tmpDivisionFlag;
        row                          = row+1;
        col                          = col-1;
    end
end

% Now that we know there is at least one free space in the Moore neighbourhood
% we assign random numbers to the empty ones
nRnd=0; sRnd=0; wRnd=0; eRnd=0; nwRnd=0; neRnd=0; seRnd=0; swRnd=0;
totalFree = checkTotalEmptyNeighbors(state,n);
n=state(row-1,col); s=state(row+1,col); w=state(row,col-1); e=state(row,col+1);
nw=state(row-1,col-1); ne=state(row-1,col+1); se=state(row+1,col+1); sw=state(row+1,col-1);
nRnd(n=='E')   = rand+(8-totalFree(row-1,col));
sRnd(s=='E')   = rand+(8-totalFree(row+1,col));
wRnd(w=='E')   = rand+(8-totalFree(row,col-1));
eRnd(e=='E')   = rand+(8-totalFree(row,col+1));
nwRnd(nw=='E') = rand+(8-totalFree(row-1,col-1));
neRnd(ne=='E') = rand+(8-totalFree(row-1,col+1));
seRnd(se=='E') = rand+(8-totalFree(row+1,col+1));
swRnd(sw=='E') = rand+(8-totalFree(row+1,col-1));

% nRnd(n=='E')   = rand+1;
% sRnd(s=='E')   = rand+1;
% wRnd(w=='E')   = rand+1;
% eRnd(e=='E')   = rand+1;
% nwRnd(nw=='E') = rand;
% neRnd(ne=='E') = rand;
% seRnd(se=='E') = rand;
% swRnd(sw=='E') = rand;

if nRnd==max([nRnd,sRnd,wRnd,eRnd,nwRnd,neRnd,seRnd,swRnd])
    state(row-1,col)=currState;
    divisionCounter(row-1,col)=currDivisionCounter;
    divisionFlag(row-1,col)=currDivisionFlag;
elseif sRnd==max([nRnd,sRnd,wRnd,eRnd,nwRnd,neRnd,seRnd,swRnd])
    state(row+1,col)=currState;
    divisionCounter(row+1,col)=currDivisionCounter;
    divisionFlag(row+1,col)=currDivisionFlag;
elseif wRnd==max([nRnd,sRnd,wRnd,eRnd,nwRnd,neRnd,seRnd,swRnd])
    state(row,col-1)=currState;
    divisionCounter(row,col-1)=currDivisionCounter;
    divisionFlag(row,col-1)=currDivisionFlag;
elseif eRnd==max([nRnd,sRnd,wRnd,eRnd,nwRnd,neRnd,seRnd,swRnd])
    state(row,col+1)=currState;
    divisionCounter(row,col+1)=currDivisionCounter;
    divisionFlag(row,col+1)=currDivisionFlag;
elseif nwRnd==max([nRnd,sRnd,wRnd,eRnd,nwRnd,neRnd,seRnd,swRnd])
    state(row-1,col-1)=currState;
    divisionCounter(row-1,col-1)=currDivisionCounter;
    divisionFlag(row-1,col-1)=currDivisionFlag;
elseif neRnd==max([nRnd,sRnd,wRnd,eRnd,nwRnd,neRnd,seRnd,swRnd])
    state(row-1,col+1)=currState;
    divisionCounter(row-1,col+1)=currDivisionCounter;
    divisionFlag(row-1,col+1)=currDivisionFlag;
elseif seRnd==max([nRnd,sRnd,wRnd,eRnd,nwRnd,neRnd,seRnd,swRnd])
    state(row+1,col+1)=currState;
    divisionCounter(row+1,col+1)=currDivisionCounter;
    divisionFlag(row+1,col+1)=currDivisionFlag;
elseif swRnd==max([nRnd,sRnd,wRnd,eRnd,nwRnd,neRnd,seRnd,swRnd])
    state(row+1,col-1)=currState;
    divisionCounter(row+1,col-1)=currDivisionCounter;
    divisionFlag(row+1,col-1)=currDivisionFlag;
end

end


