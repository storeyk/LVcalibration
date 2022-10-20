% Function to be used with 'main_CA_vary_interactions.m', which simulates
% a CA model described in Cho, Lewis, Storey, Byrne, "Designing 
% experimental conditions to use the Lotka-Volterra model to infer tumor 
% cell line interaction types", 2022.
%
% This function returns how many neighbors of a site in the CA model are 
% empty, and for proliferating cells, how many sites are occupied by cells
% of the opposite type.
%
% Author: Kathleen Storey <storeyk@lafayette.edu>
% Last revision: 10-18-2022


%%
function [totalFree,totalOpp] = checkNeighbors(state,n) 

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

%check how many of opposite type 
westOpp = zeros(n,n); eastOpp = zeros(n,n); northOpp = zeros(n,n); southOpp = zeros(n,n);
northWestOpp = zeros(n,n); northEastOpp = zeros(n,n); southWestOpp = zeros(n,n); southEastOpp = zeros(n,n);
westOpp(state=='C' & (westNghbrs=='R'|westNghbrs=='q'))=1;
westOpp(state=='R' & (westNghbrs=='C'|westNghbrs=='Q'))=1;
eastOpp(state=='C' & (eastNghbrs=='R'|eastNghbrs=='q'))=1;
eastOpp(state=='R' & (eastNghbrs=='C'|eastNghbrs=='Q'))=1;
northOpp(state=='C' & (northNghbrs=='R'|northNghbrs=='q'))=1;
northOpp(state=='R' & (northNghbrs=='C'|northNghbrs=='Q'))=1;
southOpp(state=='C' & (southNghbrs=='R'|southNghbrs=='q'))=1;
southOpp(state=='R' & (southNghbrs=='C'|southNghbrs=='Q'))=1;
northWestOpp(state=='C' & (northWestNghbrs=='R'|northWestNghbrs=='q'))=1;
northWestOpp(state=='R' & (northWestNghbrs=='C'|northWestNghbrs=='Q'))=1;
northEastOpp(state=='C' & (northEastNghbrs=='R'|northEastNghbrs=='q'))=1;
northEastOpp(state=='R' & (northEastNghbrs=='C'|northEastNghbrs=='Q'))=1;
southWestOpp(state=='C' & (southWestNghbrs=='R'|southWestNghbrs=='q'))=1;
southWestOpp(state=='R' & (southWestNghbrs=='C'|southWestNghbrs=='Q'))=1;
southEastOpp(state=='C' & (southEastNghbrs=='R'|southEastNghbrs=='q'))=1;
southEastOpp(state=='R' & (southEastNghbrs=='C'|southEastNghbrs=='Q'))=1;
totalOpp = westOpp+eastOpp+northOpp+southOpp+northEastOpp+northWestOpp...
    +southEastOpp+southWestOpp;


end

