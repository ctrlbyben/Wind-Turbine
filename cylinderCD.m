function [C_D] = cylinderCD(Re)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FUNCTION NAME: cylinderCD
%
%  PURPOSE:
%   To calculate and return the coeffecient of drag for a cyliner in cross
%   flow.  
%
%   Fit functions can be found in references:
%   
%   - White, F.M. (2006) Viscous Fluid Flow. 3rd Edition, McGraw-Hill, Boston.
%   - Nian-Sheng Cheng, *Calculation of Drag Coefficient for Arrays of ...
%      Emergent Circular Cylinders with Pseudofluid Model.*...
%      J. Hydraul. Eng. 2013.139:602-611.
%  
% INPUT
%   Re - Reynolds number [-]
%
%  OUTPUT
%   C_D - drag coeffecient [-]
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  AUTHOR: MJA
%  DATE: 2024.10.21
%
%  DESCRIPTION OF LOCAL VARIABLES
%  
%
%  FUNCTIONS CALLED
%   log10 (MATLAB)
%   tanh (MATLAB)
%
%  START OF EXECUTABLE CODE
%

% piecewise fit functions based on various reynolds number regimes for a
% cylinder in cross-flow

if Re < 2*10^5
    C_D = 11 * Re.^(-0.75) + 0.9 * (1.0 - exp(-1000./Re))...
        + 1.2 * (1.0 -exp(-(Re./4500).^0.7));

elseif Re <= 5*10^5
    C_D = 10.^(0.32*tanh(44.4504 - 8 * log10 (Re))-0.238793158);

else
    C_D = 0.1 * log10(Re) - 0.2533429;

end