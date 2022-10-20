% Function to be used with 'main_CA_vary_interactions.m', which simulates
% a CA model described in Cho, Lewis, Storey, Byrne, "Designing 
% experimental conditions to use the Lotka-Volterra model to infer tumor 
% cell line interaction types", 2022.
%
% This function returns how many neighbors of a site in the CA model are 
% empty.
%
% Author: Kathleen Storey <storeyk@lafayette.edu>
% Last revision: 10-18-2022

%%
function totalFree = checkTotalEmptyNeighbors(state,n)
% Computes the number of empty cells in the first
%order Moore neighborhood

westNghbrs = circshift(state,[0,1]);
eastNghbrs = circshift(state,[0,-1]);
northNghbrs = circshift(state,[1,0]);
southNghbrs = circshift(state,[-1,0]);
northWestNghbrs = circshift(state,[1,1]);
northEastNghbrs = circshift(state,[1,-1]);
southEastNghbrs = circshift(state,[-1,-1]);
southWestNghbrs = circshift(state,[-1,1]);

% Check how many empty cells around me
westFree = zeros(n,n); eastFree = zeros(n,n); northFree = zeros(n,n); southFree = zeros(n,n);
northWestFree = zeros(n,n); northEastFree = zeros(n,n); southWestFree = zeros(n,n); southEastFree = zeros(n,n);
westFree(westNghbrs=='E')=1;
eastFree(eastNghbrs=='E')=1;
northFree(northNghbrs=='E')=1;
southFree(southNghbrs=='E')=1;
northWestFree(northWestNghbrs=='E')=1;
northEastFree(northEastNghbrs=='E')=1;
southWestFree(southWestNghbrs=='E')=1;
southEastFree(southEastNghbrs=='E')=1;
totalFree = westFree+eastFree+northFree+southFree+northEastFree+northWestFree...
    +southEastFree+southWestFree;

end

