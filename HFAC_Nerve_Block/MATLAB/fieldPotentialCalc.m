function [V, radArray] = fieldPotentialCalc(inputCurrent, length, numCompartments, ...
    electrodeDist, electrodeLoc, sigma)
%UNTITLED2 Summary of this function goes here
%   Inputs:
%       inputCurrent: injected currrent (A)
%       

% Calculates radii
compLength = length/numCompartments;
radRight = 0:compLength:(length - (electrodeLoc * compLength));
radLeft = compLength:compLength:(length - (electrodeLoc * compLength))-compLength;
distArray = flip(radLeft);
distArray = [distArray radRight];
radArray = sqrt(electrodeDist^2 + distArray.^2);

% Calculates membrane voltage at each location
V = inputCurrent ./ (4 * pi * sigma .* radArray );
V = (V')*(1e3); %to mV

end

