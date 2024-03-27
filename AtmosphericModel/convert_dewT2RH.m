% Input:
% TD: dew point in degC
% T: air temperature in degC
% Ouput:
% RH: relative humidity in %

function RH=convert_dewT2RH(TD,T); 

%this is from TD to relative humidity (in percent)
RH=100* exp((17.625*TD)./(243.04+TD)) ./ exp((17.625*T)./(243.04+T)); 

% %this is from TD to specific humidity (q)
% %where TD is dew point in degC, p is surface presure in mb, q is specific
% %humidty in kg/kg
% e = 6.112 * exp((17.67*TD)./(TD+243.5));
% q = (0.622*e) ./(p- (0.378*e));


end