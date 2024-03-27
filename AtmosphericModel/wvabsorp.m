function z = wvabsorp(f,T,P,pv)

% function for calculating water-vapor absorption [db/km]

% f is frequency in GHz
% T is temperature in Kelvins
% P is pressure in millibars
% pv = water-vapor density in g/m^3

% Programmed 11/96 by K. Tigerstedt

%Table 5.3 in Ulably (1981)
wvline = [...
022.23515	644	1.0	2.85	1.75	0.626;
183.31012	196	41.9	2.68	2.03	0.649;
323.00000	1850	334.4	2.30	1.95	0.420;
325.15380	454	115.7	3.03	1.85	0.619;
380.19680	306	651.8	3.19	1.82	0.630;
390.00000	2199	127.0	2.11	2.03	0.330;
436.00000	1507	191.4	1.50	1.97	0.290;
438.00000	1070	697.6	1.94	2.01	0.360;
442.00000	1507	590.2	1.51	2.02	0.332;
448.00080	412	973.1	2.47	2.19	0.510];

fi = wvline(:,1);
ei = wvline(:,2);
Ai = wvline(:,3);
gi0 = wvline(:,4);
ai = wvline(:,5);
xi = wvline(:,6);

sigma = 0;

i = zeros(1,10);

for i=1:10
    %this corresponds to 5.29
    gi = gi0(i) * (P/1013) * (300/T)^xi(i) * (1+10^(-2)*ai(i)*pv*T/P);
    %this corresponds to summation term in 5.28
    sigma = sigma + Ai(i)*exp(-ei(i)/T) * gi / ( (fi(i)^2 - f^2)^2 + 4*f^2*gi^2 );
end

%this corresponds to 5.28
kh20 = 2*f^2*pv*(300/T)^2.5*sigma;

%this corresponds to 5.31
dkappa = 4.69e-6*pv*(300/T)^2.1*(P/1000)*f^2;

a1=kh20;
a2=dkappa;

%this corresponds to 5.30
z = kh20+dkappa;
return;

