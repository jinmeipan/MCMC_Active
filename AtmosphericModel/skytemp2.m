% Function for calculating the atmospheric brightness temperature
% and transmissivity
%  by K. Tigerstedt 1997
% Frequency range; typ. freqr = [1 5 10 15 19:24 30 40 50 55:65 70 80 90 100]  [GHz]
% Nadir angle; typ. thetad = 0 [deg]
% Sensor direction; up_down = 'U' for ground-based, upward looking (includes cosmic background)
%                           = 'D' for space-borne, downward looking (no cosmic background)
% step size; typ.  step = 0.05 [km]
% Ground pressure; typ. p0 = 1013 [mbar]    	
% Ground temp.; typ. T0 = 293.15 [K]
% Season (month) [1..12];
% ground-level water-vapour; typ. moist0 = 7.5 [g/m^3]
% Cloud data: cloud_data = [cloud_anal cl_lower cl_upper mv rc alfa gamma]
%       where 
%             cloud_anal  = 0/1   Inclusion of cloud analysis
%             cl_lower    = lower boundary of cloud [km]
%             cl_upper    = upper boundary of cloud [km]
%             mv          = liquid water content of cloud [g/m^3]
%             rc          = mean drop radius 
%             alfa, gamma = shape parameters

% tic

function [Tsky, tran] = skytemp2(freqr,thetad,up_down,step,p0,T0,season,moist0,cloud_data);


% % %test
% if(true)
% freqr=[6.925,10.65,18.7,36.5,89,23];
% thetad=0;
% up_down='U';
% step=0.05;
% % p0=680;
% p0=1000;
% %T0=273.15-16.8; %grouund temperature is -5 degC
% T0=273.15+20;
% season=1;
% %moist0=0.9196;
% moist0=7.5; %0.9196;
% cloud_data=[0,0,0,0,0,0,0];
% end


%# scalar i
if isempty(freqr)   % standard frequency range
   freqr = [1 5 10 15 19:24 30 40 50 55:65 70 80 90 100];
end

if isempty(cloud_data) %no cloud analysis
   cloud_data = [0 0 0 0 0 0 0];
end

cloud_anal = cloud_data(1);
cl_lower = cloud_data(2);
cl_upper = cloud_data(3);
mv = cloud_data(4);
rc = cloud_data(5);
alfa = cloud_data(6);
gamma = cloud_data(7);

ceiling = 20;  % Calculation to 20 km
Tcosmic = 2.7; % Cosmic background radiation
T0K = 273.15; % 0 deg C
theta = thetad /180 *pi;
Tsky = zeros(1,max(size(freqr)));
f = zeros(1,max(size(freqr)));
h = zeros(1,max(size(0:step:ceiling)));
i = 0;
TAU = zeros(1,max(size(freqr)));

for f = freqr
% for f=freqr(3)
   i = i+1;
   attco = 0;
   tau_z = 0;
   Tatm = 0;
   kappa_dz_old = 0;
   j=0;
  for h = ceiling:-step:0
      j=j+1;
     p = pres(p0,h);
     t = temp(T0,h,season);
     wv  = wvprof(moist0,h);
     if cloud_anal & (h >= cl_lower) & (h <= cl_upper) % Are we in the cloud?
        cloud_sca = cloudext(f*1e9, mv, rc, alfa, gamma, t-T0K, 1);  %scattering
        cloud_ext = cloudext(f*1e9, mv, rc, alfa, gamma, t-T0K, 2);  %absorption
        cloud_abs = cloud_ext -  cloud_sca;
     else
        cloud_sca = 0;
        cloud_ext = 0;
        cloud_abs = 0;
     end
     % note that we approximate t_up as t_down... see text on p.283
     % part of the integrand in 5.44
     k_o2(j)=oxabsorp(f, t, p);
     k_wv(j)=wvabsorp(f,t,p,wv);
     [z,a(j)]=wvabsorp_mod(f,t,p,wv);
     kappa_dz = ( oxabsorp(f, t, p) + wvabsorp(f,t,p,wv) ) * step / (10 * log10(exp(1))) + cloud_abs*step*1000;
     % part of the integrand in 5.44
     tau_z = tau_z + kappa_dz_old * sec(theta);  
     % this is 5.47, L_theta
     attco = exp(-tau_z);
     % this may be 5.49
     Tatm = Tatm + (1-exp(-kappa_dz*sec(theta)))*t*attco;
%      TATM(j)=Tatm;
     kappa_dz_old = kappa_dz + cloud_sca*step*1000;
  end
  total_tau_z=tau_z+kappa_dz_old*sec(theta);
  total_attco=exp(-total_tau_z);
   
  Tsky(i) = Tatm;
  tran(i) = total_attco;
  TAU(i) = total_tau_z;
end

end
