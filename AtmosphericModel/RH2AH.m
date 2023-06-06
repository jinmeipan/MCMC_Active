function ah=RH2AH(rh,Tair);

%calcualte absolute humidity (g/m^3) from temperature (Tair, in K) and
%relative humidity (in %)

% inputs:
% Tair, air temperature in K. Exmple: Tair=273+20;
% rh, relative humidity in %. Example: rh=80;

%Saturated Vapour pressure in Pa, Pws
%      A         m         Tn   max-error  Temperature-range 
data=[6.116441 7.591386 240.7263 0.083% -20...+50°C
    6.004918 7.337936 229.3975 0.017% +50...+100°C
    5.856548 7.27731 225.1033 0.003% +100...+150°C
    6.002859 7.290361 227.1704 0.007% +150...+200°C
    9.980622 7.388931 263.1239 0.395% +200...+350°C
    6.089613 7.33502 230.3921 0.368% 0...+200°C
    6.114742 9.778707 273.1466 0.052];% -70...0°C

row=1;   %pick the first row in the Look-Up table
A=data(row,1);
m=data(row,2);
Tn=data(row,3);

T=Tair-273.15;
Pws = A*10.^(m.*T ./ (T+Tn));
Pws = Pws*100;  %convert from hPa to Pa

%Vapour pressure in Pa, Pw
Pw = Pws .* rh/100;


%aboslute humidity
C= 2.16679; %gK/J
ah = C * Pw ./ Tair;


end