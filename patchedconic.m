function [tot_time, dv_tot] = patchedconic (Bennu,Earth,Launch,Survey,Approach,theta_launch,u_s,u_e,u_b)
%Alex Philpott

%% Calculate Solar Frame Departure and Arrival dvs

%Hohmann Transfer burn calculation
ra = Earth.a; %Depart at average Earth radius
rp = Bennu.P;
at = (rp + ra)./2;
et = (ra-rp)./(ra+rp);

vde = sqrt(u_s.*(2./ra - 1./Earth.a));
vdt = sqrt(u_s.*(2./ra - 1./at));
dvd = abs(sqrt(vdt.^2 + vde.^2 - 2*vdt*vde*cosd(abs(Earth.i-Bennu.i)))); %Transfer insertion and inclination change
dva = abs(sqrt(u_s.*(2./rp - 1./at)) - sqrt(u_s.*(2./rp - 1./Bennu.a)));

figure; hold on; axis equal;
set(gca,'DefaultLineLineWidth',2);
xlabel('(km)');
ylabel('(km)');
theta = linspace(0,2*pi);
theta_t = linspace(pi,2*pi);
rf_E = Earth.a;
rf_B = Bennu.a.*(1-Bennu.e.^2)./(1+Bennu.e.*cos(theta));
rf_t = at.*(1-et.^2)./(1+et.*cos(theta_t));
plot(rf_E.*cos(theta),rf_E.*sin(theta),'b','DisplayName','Earth');
plot(rf_B.*cos(theta),rf_B.*sin(theta),'k','DisplayName','Bennu');
plot(rf_t.*cos(theta_t),rf_t.*sin(theta_t),'g','DisplayName','Departure Hohmann Transfer');

%% Calculate departure burn
v_infin1 = dvd;
vp = sqrt(2.*u_e./Launch.P+ v_infin1.^2);
dv1 = vp-Launch.vp;
Transfer_h.h = Launch.P.*vp;
Transfer_h.vp = vp;
Transfer_h.P = Launch.P;
Transfer_h.e = sqrt((v_infin1.*Transfer_h.h./u_e).^2 + 1);
Transfer_h.thetainfin = acos(-1./Transfer_h.e);
Transfer_h.beta = pi - Transfer_h.thetainfin;

fprintf('Depart launch orbit with delta-v %.3f km/s. Burn at perigee.\n',dv1);

%% Calculate arrival burn
v_infin2 = dva;
vp = sqrt(2.*u_b./Survey.r + v_infin2.^2);
dv2 = Survey.v - vp;
Arrival_h.h = Survey.r.*vp;
Arrival_h.vp = vp;
Arrival_h.P = Survey.r;
Arrival_h.e = sqrt((v_infin2.*Arrival_h.h./u_b).^2+1);
Arrival_h.thetainfin = acos(-1./Arrival_h.e);
Arrival_h.beta = pi - Arrival_h.thetainfin;

fprintf('Establish orbit around Bennu with delta-v %.3f km/s. Burn at periapsis.\n',dv2);

%% Calculate departure time
%Finding TOF
T = 2*pi*sqrt(at.^3/u_s);
Tda = 0.5.*T; %sec
%fprintf('The journey to Bennu will take %.3f days.\n',Tda./24./3600);

%Finding departure angle
phi_d = pi - Bennu.n.*0.5.*T %This is the true anomaly of Bennu when the object should be launched from Earth

%Calculate position of Bennu at launch time
dt = Launch.JD - Bennu.tp; %days
dt = dt.*24.*3600; %sec
M = Bennu.n.*dt;
fh = @(x) x - Bennu.e.*sin(x) - M;
f_prime = @(x) 1 - Bennu.e.*cos(x);

E = newton_method(fh,f_prime,pi/4,1e-6); %rad
E_frac = E./2/pi;
E_frac = E_frac - floor(E_frac);
E_frac = E_frac.*2.*pi; %this is the position of Bennu in its current orbit at launch time

theta_b = 2*atan(sqrt((1+Bennu.e)/(1-Bennu.e))*tan(E_frac/2));
theta_e = theta_b - theta_launch + 2*pi; %defining position of Earth based on position of Bennu

dt = 0.5.*24*3600; %0.5 day in sec
[tolerances,~,~] = sim_angle(Earth,Bennu,f_prime,theta_e,theta_b,dt,phi_d,0);

[minimum,I] = min(tolerances);

waiting_time = dt.*I;
fprintf('Wait %.3f days from Sep. 8 2016 for Earth and Bennu to be aligned, then launch. Arrive at Bennu %.3f days later.\n',waiting_time/24/3600,...
    Tda/24/3600);
fprintf('Arrival on Bennu occurs %.3f days later.\n',(Tda+waiting_time)/24/3600);

deltat = Tda;
[~,theta_e,theta_b] = sim_angle(Earth,Bennu,f_prime,theta_e,theta_b,dt,0,deltat);

fprintf('At arrival, the true anomaly of Bennu is %.3f rad and the true anomaly of Earth is %.3f.\n',theta_b,theta_e);

%% Calculating Surveying Time
phi_ret = Earth.n.*Tda - pi
[tolerances,~,~] = sim_angle(Earth,Bennu,f_prime,theta_e,theta_b,dt,phi_ret);

[minimum,I] = min(tolerances);

return_waiting_time = dt.*I;
fprintf('Depart from Bennu %.3f days later and arrive at Earth %.3f days after departure.\n',...
    return_waiting_time/3600/24,Tda/3600/24);
[~,theta_e,theta_b] = sim_angle(Earth,Bennu,f_prime,theta_e,theta_b,dt,0,return_waiting_time);
fprintf('At departure from Bennu, the true anomaly of Bennu is %.3f and the true anomaly of Earth is %.3f.\n',...
    theta_b,theta_e);

%% Calculating return delta vs in Solar Frame
ra = Bennu.Q;
rp = Earth.a;
at = (ra + rp)/2;
et = (ra-rp)./(ra+rp);
theta_t = linspace(pi,2*pi);
rft2 = at.*(1-et.^2)./(1+et.*cos(theta_t));
plot(rft2.*cos(theta_t),rft2.*sin(theta_t),'r','DisplayName','Return Hohmann Transfer');
legend('Location','Northwest');
print('patchedconic','-dpng');

vdb = sqrt(u_s.*(2./ra - 1./Bennu.a));
vdt = sqrt(u_s.*(2./ra - 1./at));
dvd = abs(sqrt(vdt.^2 + vdb.^2 - 2*vdt*vdb*cosd(abs(Earth.i-Bennu.i)))); %Transfer insertion and inclination change
dva = abs(sqrt(u_s.*(2./rp - 1./at)) - sqrt(u_s./rp));



%% Calculating Bennu departure delta v
v_infin = dvd;
vp = sqrt(2.*u_b./Survey.r+ v_infin.^2);
dv3 = vp - Survey.v;

v_infin2 = dva;
vp = sqrt(2.*u_e./Approach.r + v_infin2.^2);
dv4 = 0;

fprintf('Depart Bennu orbit with delta-v %.3f.\n',dv3);

tot_time = waiting_time + T + return_waiting_time;
dv_tot = abs(dv1) + abs(dv2) + abs(dv3) + abs(dv4);
end