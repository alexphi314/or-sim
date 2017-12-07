%Alex Philpott Final Project
clc;
clear;
close all;
mode = 2;

%% Defining Constants
G = 6.67e-11./1000^3; %km3/kg/s2
u_e = G.*5.974e24; %km3/s2
u_s = G.*1.989e30; %km3/s2
u_b = G.*6e10; %km3/s2
RE = 6371.008; %km

Sun_mass = 1.989e30; %kg
OR_mass = 2110; %kg

au = 149597870700./1000; %km

%Defining Bennu orbital parameters
Bennu.e = .2037451095574896;
Bennu.a = 1.126391026024589.*au; %km
Bennu.i = 6.034939391321328; %deg
Bennu.RAAN = 2.060867570687039; %deg
Bennu.wp = 66.22306847249314; %deg
Bennu.M = 101.7039478618686; %deg
Bennu.tp = 2455439.141945873525; %JD
Bennu.T = 436.6487281874828.*24.*3600; %sec
Bennu.n = .8244613501895456.*pi./(180.*24.*3600); %rad/s
Bennu.Q = 1.355887689026542.*au; %km %r_a
Bennu.P = -Bennu.Q*(Bennu.e-1)/(Bennu.e+1); %km
Bennu.mass = 6e10; %kg

%Defining Earth orbital parameters
Earth.e = 0.01671123;
Earth.a = 1.00000261.*au; %km
Earth.i = -0.00001531; %deg
Earth.RAAN = 0; %deg
Earth.wp = 102.93768193; %deg
Earth.T = 365.256.*24.*3600; %sec
Earth.n = sqrt(u_s./Earth.a.^3); %rad/s
Earth.P = 147.09e6; %km
Earth.Q = 152.10e6; %km
Earth.mass = 5.97219e24; %kg

%Defining launch orbit parameters
Launch.P = 113.09.*1.60934 + RE; %km
Launch.Q = -3443.92.*1.60934 - RE; %km
Launch.a = (Launch.P + Launch.Q)./2; %km
Launch.i = 29.447; %deg
Launch.e = (Launch.Q-Launch.P)./(Launch.Q + Launch.P);
Launch.h = sqrt(u_e.*Launch.P.*(1+Launch.e));
Launch.vp = Launch.h./Launch.P;
Launch.JD = 2457639.5; %JD

%Defining Jupiter orbit parameters
Jupiter.a = 5.20288700.*au; %km
Jupiter.e = 0.04838624;
Jupiter.i = 1.30439695; %deg
Jupiter.wp = 14.72847983; %deg
Jupiter.RAAN = 100.47390909; %deg
Jupiter.radius = 69911; %km
Jupiter.mass = 1.8981e27; %kg
Jupiter.T = 11.86.*365.25.*24.*3600; %sec

%Defining survey orbit parameters
Survey.r = 5; %km
Survey.i = 90; %deg
Survey.v = sqrt(u_b./Survey.r);

%Defining final fly-by orbit parameters
Approach.r = 200;

%Calculating angle between Earth and Bennu at launch
earth_r = [1.459452879414980E+08,-3.761408174440899E+07,2.573073861077428E+02]; %km
bennu_r = [-2.479434376975826E+06,1.351546412664345E+08,1.428915405096796E+07]; %km
theta_launch = acos(dot(earth_r,bennu_r)/norm(earth_r)/norm(bennu_r));

%% 2-Body Patched Conic Total delta-V calculation
[tot_time, dv_tot] = patchedconic(Bennu,Earth,Launch,Survey,Approach,theta_launch,u_s,u_e,u_b);
fprintf('Total time from launch to return: %.3f days. Total dv: %.3f km/s.\n',tot_time/24/3600,dv_tot);

%% 4-Body (Jupiter, Earth, Bennu, Sun) delta-V calculation
fprintf('\n4 Body Simulations \n');
phasing_duration = 379.*24.*3600; %sec
[Sun_pos,Earth_pos,Jupiter_pos,Bennu_pos,OR_pos,t,OR_veloc] = orbit_sim(Bennu,Earth,Launch,Survey,Approach,phasing_duration,u_s,u_e,u_b);

%% Parsing or Loading Ephemerides
if mode == 1 %parse
    [E_Sun,VE_Sun,E_Earth,VE_Earth,E_Jupiter,VE_Jupiter,E_Bennu,VE_Bennu,E_OR,VE_OR] = parse_jpl();
    %saving variables
    save('parsed_files','E_Sun','VE_Sun','E_Earth','VE_Earth','E_Jupiter','VE_Jupiter','E_Bennu','VE_Bennu','E_OR','VE_OR');
else
    load('parsed_files.mat');
end
%% Plotting
Sun_x0 = E_Sun(1,:); %km
Sun_v0 = VE_Sun(1,:); %km/s

Earth_x0 = E_Earth(1,:); %km
Earth_v0 = VE_Earth(1,:); %km/s

Jupiter_x0 = E_Jupiter(1,:); %km
Jupiter_v0 = VE_Jupiter(1,:); %km/s

Bennu_x0 = E_Bennu(1,:); %km
Bennu_v0 = VE_Bennu(1,:); %km/s

OR_x0 = E_OR(1,:); %km
OR_v0 = VE_OR(1,:); %km/s

hold on; axis equal;
%h1 = plot(Sun_pos(:,1),Sun_pos(:,2),'r','DisplayName','Sun');
%plot(Sun_x0(1),Sun_x0(2),'r*');
%plot(Sun_pos(end,1),Sun_pos(end,2),'ro');
h2 = plot(Earth_pos(:,1),Earth_pos(:,2),'b','DisplayName','Earth');
h3 = plot(E_Earth(1:379,1),E_Earth(1:379,2),':r','DisplayName','JPL Earth');
%plot(Earth_x0(1),Earth_x0(2),'b*');
plot(Earth_pos(end,1),Earth_pos(end,2),'r*');
%plot(Jupiter_pos(:,1),Jupiter_pos(:,2),'r');
%plot(E_Jupiter(1:379,1),E_Jupiter(1:379,2),'--b');
%h3 = plot(Bennu_pos(:,1),Bennu_pos(:,2),'k','DisplayName','Bennu');
%plot(E_Bennu(1:379,1),E_Bennu(1:379,2),'--r');
%plot(Bennu_x0(1),Bennu_x0(2),'k*');
%plot(Bennu_pos(end,1),Bennu_pos(end,2),'ko');
h4 = plot(OR_pos(:,1),OR_pos(:,2),'g','DisplayName','Osiris Rex');
h5 = plot(E_OR(1:379,1),E_OR(1:379,2),'--k','DisplayName','JPL Osiris Rex');
plot(OR_x0(1),OR_x0(2),'go');
plot(OR_pos(end,1),OR_pos(end,2),'g*');
legend([h2 h3 h4 h5],'Location','Northwest');

print('orbits','-dpng');

%% Analyzing
r = zeros(size(OR_pos,1),1);
for k = 1:size(OR_pos,1)
    r(k) = norm(OR_pos(k) - Earth_pos(k));
end
test = r < 1e4;
v = [];
for k = 1:size(OR_veloc,1)
    v(k,1) = norm(VE_OR(k,:));
    v(k,2) = norm(OR_veloc(k,:));
end