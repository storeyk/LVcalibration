
% This function contains the model used in Cho,
% Lewis, Storey, Byrne, "Designing experimental conditions to use the
% Lotka-Volterra model to infer tumor cell line interaction types", 2022.
%
% Author: Allison Lewis <lewisall@lafayette.edu>
% Last revision: 10-18-2022

function dval = tumorTwoComp( t, val, Params ) 


rS = Params(1);
rR = Params(2);
KS = Params(3);
KR = Params(4);
gammaS = Params(5); 
gammaR = Params(6); 

S = val(1);
R = val(2);


dval(1,1) =  rS * S * ( 1 - S/KS - gammaR*R/KS);
dval(2,1) =  rR * R * ( 1 - R/KR - gammaS*S/KR);


end 


