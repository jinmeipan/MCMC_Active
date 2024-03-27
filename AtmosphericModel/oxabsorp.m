function z = oxabsorp(f,T,P)

% Function for calculating oxygen absorption coefficient [dB/km]

% f is frequency in GHz
% T is temperature in Kelvins
% P is pressure in millibars

% programmed 11/96 by K. Tigerstedt

% load oxfic.par
% Absorption data now included in function 02/98

%Table 5.4 in Ulably (1981)
oxfic = [...
56.2648	118.7503	4.51e-4	-2.14e-5;
58.4466	62.4863	4.94e-4	-3.78e-4;
59.5910	60.3061	3.52e-4	-3.92e-4;
60.4348	59.1642	1.86e-4	-2.68e-4;
61.1506	58.3239	3.30e-5	-1.13e-4;
61.8002	57.6125	-1.03e-4	3.44e-5;
62.4112	56.9682	-2.23e-4	1.65e-4;
62.9980	56.3634	-3.32e-4	2.84e-4;
63.5685	55.7838	-4.32e-4	3.91e-4;
64.1278	55.2214	-5.26e-4	4.93e-4;
64.6789	54.6711	-6.13e-4	5.84e-4;
65.2241	54.1300	-6.99e-4	6.76e-4;
65.7647	53.5957	-7.74e-4	7.55e-4;
66.3020	53.0668	-8.61e-4	8.47e-4;
66.8367	52.5422	-9.11e-4	9.01e-4;
67.3694	52.0212	-1.03e-3	1.03e-3;
67.9007	51.5030	-9.87e-4	9.86e-4;
68.4308	50.9873	-1.32e-3	1.33e-3;
68.9601	50.4736	-7.07e-4	7.01e-4;
69.4887	49.9618	-2.58e-3	2.64e-3];

fnp = oxfic(:,1);
fnm = oxfic(:,2);
Ynp = oxfic(:,3);
Ynm = oxfic(:,4);

sigma = 0;

N = zeros(1,max(size(1:2:39)));

for N=1:2:39

 i = (N-1)/2+1;  % index for parameters from file oxfic.par

%these four are from 5.37-5.39 
 gaN = 1.18*(P/1013)*(300/T)^0.85;
 gab = 0.49*(P/1013)*(300/T)^0.89;
 dnp = sqrt(N*(2*N+3)/((N+1)*(2*N+1)));
 dnm = sqrt((N+1)*(2*N-1)/(N*(2*N+1)));

 %these correspond to 5.35
 gNplusf = (gaN*dnp^2 + P*(f-fnp(i))*Ynp(i)) / ( (f-fnp(i))^2 + gaN^2 );
 gNplusmf = (gaN*dnp^2 + P*(-f-fnp(i))*Ynp(i)) / ( (-f-fnp(i))^2 + gaN^2 );
 gNminusf = (gaN*dnm^2 + P*(f-fnm(i))*Ynm(i)) / ( (f-fnm(i))^2 + gaN^2 );
 gNminusmf = (gaN*dnm^2 + P*(-f-fnm(i))*Ynm(i)) / ( (-f-fnm(i))^2 + gaN^2 );

 %this corresponds to 5.36
 thetaN = 4.6e-3*(300/T)*(2*N+1)*exp(-6.89e-3*N*(N+1)*(300/T));

 % this corresponds to the second term on RHS of 5.34
 sigma = sigma + thetaN*(gNplusf + gNplusmf + gNminusf + gNminusmf);

end
%this corresponds to 5.34
Fdot = 0.7*gab/(f^2 + gab^2) + sigma;

%this corresponds to 5.33
z = 1.61e-2*f^2*(P/1013)*(300/T)^2*Fdot;

return;
