function z = pres(p0,h)

% Function for calculating vertical pressure profiles
% Assume exponetial pressure profile
% p0 = sea level pressure
% h = height [km]

Hp = 7.7 ; % Pressure scale height; 7.7 km

z = p0*exp(-h / Hp);

return;
