% Function to be used with 'main_CA_vary_interactions.m', which simulates
% a CA model described in Cho, Lewis, Storey, Byrne, "Designing 
% experimental conditions to use the Lotka-Volterra model to infer tumor 
% cell line interaction types", 2022.
%
% This function computes the volume of a CA spheroid. It assumes that the
% spheroid is an ellipsoid, with the formula vol=0.5*hor*ver^2, where
% hor is the length in the horizontal direction through the spheroid center
% and ver is the length in the vertical direction.
%
% Author: Kathleen Storey <storeyk@lafayette.edu>
% Last revision: 10-18-2022


%%
function vol = computeVolume(state,n)

CELL_SIZE = 0.018;

vol=zeros(1,9);


% Compute distances in horizontal and vertical directions through the
% centre
tmp                    = zeros(1,n);
tmp(state(:,:)=='C') = 1;
numP = sum(tmp);
tmp                    = zeros(1,n);
tmp(state(:,:)=='R') = 1;
numR = sum(tmp);
tmp                    = zeros(1,n);
tmp(state(:,:)=='Q') = 1;
numQ = sum(tmp);
tmp                    = zeros(1,n);
tmp(state(:,:)=='q') = 1;
numq = sum(tmp);
tmp                    = zeros(1,n);
tmp(state(:,:)=='N'|state(:,:)=='I') = 1;
numN = sum(tmp);
tmp                    = zeros(1,n);
tmp(state(:,:)=='n') = 1;
numn = sum(tmp);

total=sum([numP, numR, numQ, numq, numN, numn]);

vol(3:9)=[numP numR numQ numq numN numn total];

% Compute necrotic volume
tmp                    = zeros(1,n);
tmp(state(n/2,:)=='N'|state(n/2,:)=='n'|state(n/2,:)=='I') = 1;
%hor_nec=sum(tmp)
hor                    = sum(tmp)*CELL_SIZE;

tmp                    = zeros(n,1);
tmp(state(:,n/2)=='N'|state(:,n/2)=='n'|state(:,n/2)=='I') = 1;
ver                    = sum(tmp)*CELL_SIZE;

nec_vol = 0.5*hor*ver^2;

% Compute total volume
tmp                    = zeros(1,n);
tmp(state(n/2,:)~='E') = 1;
hor                    = sum(tmp)*CELL_SIZE;

tmp                    = zeros(n,1);
tmp(state(:,n/2)~='E') = 1;
ver                    = sum(tmp)*CELL_SIZE;

%total volume
total_vol=0.5*hor*ver^2;
vol(1) = total_vol;

%proportion of necrotic volume
if total_vol~=0
    vol(2) = nec_vol/total_vol;
else
    vol(2) = 0;
end



