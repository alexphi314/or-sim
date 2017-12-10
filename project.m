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

%Defining timeline of misson;
ts = [0; %Sep 9 2016 // Launch
    111; %Dec 29 2016 // Course Adjustment
    378; %Sep 22 2017 // Earth Fly-by
    540; %Simulation Adjustment
    754; % Course Adjustment
    768; %Course Adjustment
    800; %Aug 17 2018 // Bennu Arrival %Originally 707
    1636; %Mar 3 2021 // Bennu Departure
    2571; %Sep 24 2023 // Earth Fly-by
    ];
    
%% 2-Body Patched Conic Total delta-V calculation
[tot_time, dv_tot] = patchedconic(Bennu,Earth,Launch,Survey,Approach,theta_launch,u_s,u_e,u_b);
fprintf('Total time from launch to return: %.3f days. Total dv: %.3f km/s.\n',tot_time/24/3600,dv_tot);

%% 4-Body (Jupiter, Earth, Bennu, Sun) delta-V calculation
fprintf('\n4 Body Simulations \n');
phasing_duration = 379.*24.*3600; %sec
[Sun_pos,Earth_pos,Jupiter_pos,Bennu_pos,OR_pos,t,OR_veloc] = orbit_sim(Bennu,Earth,Launch,Survey,Approach,phasing_duration,u_s,u_e,u_b,ts);

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

figure;axis equal;  hold on; %view(3);
plot3(Earth_pos(:,1),Earth_pos(:,2),Earth_pos(:,3),'b','DisplayName','Earth');
h4 = plot3(OR_pos(1:ts(3),1),OR_pos(1:ts(3),2),OR_pos(1:ts(3),3),'g','DisplayName','Osiris Rex');
h5 = plot3(E_OR(1:ts(3),1),E_OR(1:ts(3),2),E_OR(1:ts(3),3),'--k','DisplayName','JPL Osiris Rex');
legend('Location','Northwest');

figure; hold on; axis equal;
plot3(OR_pos(ts(3):end,1),OR_pos(ts(3):end,2),OR_pos(ts(3):end,3),'g','DisplayName','Osiris Rex');
plot3(E_OR(ts(3):ts(7),1),E_OR(ts(3):ts(7),2),E_OR(ts(3):ts(7),3),'--k','DisplayName','JPL Osiris Rex');
plot3(OR_pos(602,1),OR_pos(602,2),OR_pos(602,3),'*g');
plot3(E_OR(602,1),E_OR(602,2),E_OR(602,3),'*k');
plot3(OR_pos(754,1),OR_pos(754,2),OR_pos(754,3),'*g');
plot3(OR_pos(end,1),OR_pos(end,2),OR_pos(end,3),'go');
plot3(E_OR(800,1),E_OR(800,2),E_OR(800,3),'ko');
%plot3(Bennu_pos(ts(3):end,1),Bennu_pos(ts(3):end,2),Bennu_pos(ts(3):end,3),'r','DisplayName','Bennu');
plot3(Bennu_pos(end,1),Bennu_pos(end,2),Bennu_pos(end,3),'ro');
plot3(E_Bennu(800,1),E_Bennu(800,2),E_Bennu(800,3),'k*');
%plot3(E_Bennu(ts(3):ts(7),1),E_Bennu(ts(3):ts(7),2),E_Bennu(ts(3):ts(7),3),'--r','DisplayName','JPL Bennu');

% figure; hold on; axis equal;
% plot(E_OR(:,1),E_OR(:,2),'g');
% plot(E_Earth(:,1),E_Earth(:,2),'b');
% plot(E_Bennu(:,1),E_Bennu(:,2),'k');

%% Analyzing
r = zeros(size(OR_pos,1),1);
for k = 1:size(OR_pos,1)
    r(k,1) = norm(OR_pos(k,:));
    r(k,2) = norm(E_OR(k,:));
    r(k,3) = abs((r(k,1)-r(k,2))/r(k,2).*100);
end
v = [];
for k = 1:size(OR_veloc,1)
    v(k,1) = norm(VE_OR(k,:));
    v(k,2) = norm(OR_veloc(k,:));
    v(k,3) = abs((v(k,1)-v(k,2))/v(k,1).*100);
    if k > 1
        v(k,4) = abs(v(k,1)-v(k-1,1))/v(k-1,1).*100;
        if v(k,4) > 0.4
            fprintf('delta-V %% of %.3f found on Day %i\n',v(k,4),k);
        end
    else
        v(k,4) = 0;
    end
end

%Grabbing vectors
before_burn = VE_OR(753,:);
after_burn = VE_OR(754,:);
bb_u = before_burn./norm(before_burn);
ab_u = after_burn./norm(after_burn);

test = t./24./3600;

r_arrival = [];
for k = 1:size(OR_pos,1)
    r_b = OR_pos(k,:) - E_Bennu(k,:);
    r_arrival(k,1) = norm(r_b);
    r_arrival(k,2) = norm(E_OR(k,:) - E_Bennu(k,:));
end
[min_dist tca] = min(r_arrival(:,1));
[min_dist_jpl tca2] = min(r_arrival(:,2));
fprintf('On Day %i Osiris Rex is %.3f km from Bennu.\n',ts(4),norm(r_arrival));

bennu_error = [];
for k = 1:size(Bennu_pos,1)
    bennu_error(k,1) = norm(Bennu_pos(k,:));
    bennu_error(k,2) = norm(E_Bennu(k,:));
    bennu_error(k,3) = abs(bennu_error(k,2)-bennu_error(k,1))/bennu_error(k,2)*100;
end

movement_u = OR_pos(540,:) - Bennu_pos(540,:);
movement_u = movement_u./norm(movement_u);