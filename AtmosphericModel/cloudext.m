% Function for calculating volume scattering/extinction coefficients 
% for hydrometeors using Mie scattering efficiencies
% Programmed by K.Tigerstedt 

function K = cloudext(f, mv, rc, alfa, gam, t, typ)

% f = frequency
% mv = mass density of liquid water
% rc = mode radius of drop-size distribution
% alfa,gamma = shape parameters
% t = temperature (deg C)
% typ = type : 1 = scattering, 2 = extinction



if t >= 0
   epsr = epsws(f,0,t);    % Dielectric coefficient of water
else
   epsr = eps_ice(f,t);    % Dielectric coefficient of ice
end

c = 3e8;

lr = 1e-7;     % Lower boundary of drop radii (m)
ur = 1e-3;      % Upper boundary of drop radii (m)
trace = [];
tol = 10e-3;

lX = 2*pi*lr*f/c;
uX = 2*pi*ur*f/c;

K = (c^3/(8*pi^2*f^3)) * quadg_m('cloudint',lX,uX,tol,trace,f,mv,rc,alfa,gam,epsr,typ);

return;

