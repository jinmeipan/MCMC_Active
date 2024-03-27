function z = temp(T0,h,season)

% Function for calculating vertical temperature profile
% T0 = ground temperature [K]
% h = height [km]
% season = 'summer' or 'winter' or month [1..12]
% programmed 3/97 by K.Tigerstedt
% updated 12/97

h1 = 8;

if season == 'winter'
        k1 = 2.96e-3;
        k2 = -5.83e-3;
        k3 = 3.77e-4;
        T2 = 217;
        h2 = 8.5;
elseif season == 'summer'
        k1 = -3.21e-2;
        k2 = 3.69e-3;
        k3 = -3.32e-4;
        T2 = 225;
        h2 = 10.0;
elseif season == 1
        k1 = 2.958e-3;
        k2 = -5.832e-3;
        k3 = 3.771e-4;
        T2 = 217;
        h2 = 8.5;
elseif season == 2
        k1 = -2.196e-3;
        k2 = -4.419e-3;
        k3 = 2.673e-4;
        T2 = 218.33;
        h2 = 8.75;
elseif season == 3
        k1 = -1.717e-2;
        k2 = -4.056e-4;
        k3 = -4.252e-5;
        T2 = 219.67;
        h2 = 9;
elseif season == 4
        k1 = -2.702e-2;
        k2 = 2.327e-3;
        k3 = -2.521e-4;
        T2 = 221.00;
        h2 = 9.25;
elseif season == 5
        k1 = -3.33e-2;
        k2 = 3.915e-3;
        k3 = -3.616e-4;
        T2 = 222.33;
        h2 = 9.5;
elseif season == 6
        k1 = -3.461e-2;
        k2 = 4.154e-3;
        k3 = -3.632e-4;
        T2 = 223.67;
        h2 = 9.75;
elseif season == 7
        k1 = -3.213e-2;
        k2 = 3.69e-3;
        k3 = -3.317e-4;
        T2 = 225;
        h2 = 10.0;
elseif season == 8
        k1 = -3.05e-2;
        k2 = 3.423e-3;
        k3 = -3.182e-4;
        T2 = 223.67;
        h2 = 9.75;
elseif season == 9
        k1 = -2.928e-2;
        k2 = 3.543e-3;
        k3 = -3.472e-4;
        T2 = 222.33;
        h2 = 9.5;
elseif season == 10
        k1 = -1.997e-2;
        k2 = 1.777e-3;
        k3 = -2.621e-4;
        T2 = 221;
        h2 = 9.25;
elseif season == 11
        k1 = -1.533e-2;
        k2 = 8.062e-4;
        k3 = -2.24e-4;
        T2 = 219.67;
        h2 = 9.0;
elseif season == 12
        k1 = -5.856e-3;
        k2 = -2.773e-3;
        k3 = 1.086e-4;
        T2 = 218.33;
        h2 = 8.75;
else
        disp('Error in seasonal input');
end


T1 = T0*(1 + k1*h1 + k2*h1^2 + k3*h1^3);

if h < 8
        z = T0*(1 + k1*h + k2*h^2 + k3*h^3);
elseif (h >= 8) & (h < h2)
        z = (T1-T2)*(h-h2) / (h1 - h2) + T2;
elseif (h >= h2) & (h <= 20)
        z = T2;
end

return;
