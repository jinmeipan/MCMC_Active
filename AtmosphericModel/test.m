% Absolute humidity in the atmosphere ranges from near zero to roughly 30 grams per cubic meter when the air is saturated at 30 °C.[4]
% moist0=1:0.5:30;
% rh=20:20:90;
rh=90;
Tair=273.15-10;
moist=RH2AH(rh,Tair);
    
tb=nan(length(moist),6);

for i=1:length(rh)
    freqr=[6.925,10.65,18.7,36.5,89,23];
    thetad=55;
    up_down='U';
    step=0.05;
    p0=1013;
    T0=273.15-5; %grouund temperature is -5 degC
    season=1;
    cloud_data=[0,0,0,0,0,0,0];
    
    [Tsky, tran] = skytemp2(freqr,thetad,up_down,step,p0,T0,season,moist(i),cloud_data);
    tb(i,:)=Tsky;
end

tb

%%
figure;
plot(rh,tb(:,1));hold on;
plot(rh,tb(:,2));
plot(rh,tb(:,3));
plot(rh,tb(:,4));
plot(rh,tb(:,5));

figure;
plot(rh,tb(:,6)-tb(:,5),'kx-');