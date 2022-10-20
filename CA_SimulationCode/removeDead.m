% Function to be used with 'main_CA_vary_interactions.m', which simulates
% a CA model described in Cho, Lewis, Storey, Byrne, "Designing 
% experimental conditions to use the Lotka-Volterra model to infer tumor 
% cell line interaction types", 2022.
%
% This function removes a dead cell from the CA simulation and shifts 
% a cell chain towards outermost cell with largerst distance to the center.
%
% Author: Kathleen Storey <storeyk@lafayette.edu>
% Last revision: 10-18-2022

%%
function [state,divisionCounter,removeFlag] = removeDead(state,divisionCounter,removeFlag,row,col,n)

% Find non-empty cell with largest distance to the center
h=n/2;
rSquared = zeros(n);
for iRow=1:n
    for iCol=1:n
        if state(iRow,iCol)~='E'
            rSquared(iRow,iCol)=(iRow-h)^2+(iCol-h)^2;
        end
    end
end

[wantedRow,wantedCol] = find(rSquared==max(max(rSquared)));
% In case there are more than one cell with the same distance take the
% first one encountered
wRow = wantedRow(1);
wCol = wantedCol(1);


% Shift from (row,col) to (wRow,wCol)
while row~=wRow || col~=wCol
    if row>wRow && col>wCol
        state(row,col) = state(row-1,col-1);
        divisionCounter(row,col) = divisionCounter(row-1,col-1);
        removeFlag(row,col) = removeFlag(row-1,col-1);
        row = row-1; col = col-1;
    elseif row<wRow && col<wCol
        state(row,col) = state(row+1,col+1);
        divisionCounter(row,col) = divisionCounter(row+1,col+1);
        removeFlag(row,col) = removeFlag(row+1,col+1);
        row = row+1; col = col+1;
    elseif row>wRow && col<wCol
        state(row,col) = state(row-1,col+1);
        divisionCounter(row,col) = divisionCounter(row-1,col+1);
        removeFlag(row,col) = removeFlag(row-1,col+1);
        row = row-1; col = col+1;
    elseif row<wRow && col>wCol
        state(row,col) = state(row+1,col-1);
        divisionCounter(row,col) = divisionCounter(row+1,col-1);
        removeFlag(row,col) = removeFlag(row+1,col-1);
        row = row+1; col = col-1;
    elseif row<wRow && col==wCol
        state(row,col) = state(row+1,col);
        divisionCounter(row,col) = divisionCounter(row+1,col);
        removeFlag(row,col) = removeFlag(row+1,col);
        row = row+1;
    elseif row>wRow && col==wCol
        state(row,col) = state(row-1,col);
        divisionCounter(row,col) = divisionCounter(row-1,col);
        removeFlag(row,col) = removeFlag(row-1,col);
        row = row-1;
    elseif row==wRow && col<wCol
        state(row,col) = state(row,col+1);
        divisionCounter(row,col) = divisionCounter(row,col+1);
        removeFlag(row,col) = removeFlag(row,col+1);
        col = col+1;
    elseif row==wRow && col>wCol
        state(row,col) = state(row,col-1);
        divisionCounter(row,col) = divisionCounter(row,col-1);
        removeFlag(row,col) = removeFlag(row,col-1);
        col = col-1;
    end
end
state(row,col) = 'E';
divisionCounter(row,col) = 0;
removeFlag(row,col) = 0;
    
end