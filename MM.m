function [output] = MM(conc,Vmax,Km)

% function to calculate Michaelis-Menten equation at given: 
% conc = substrate concentration
% Vmax = maximum uptake velocity per min
% Km = michaelis-menten constant in [g/L] == concentration at half saturation

output = Vmax*conc/(Km+conc) ; 

end 