function [Sun_pos,Earth_pos,Jupiter_pos,Bennu_pos,OR_pos,t,OR_vel] = orbit_sim(Bennu,Earth,Launch,Survey,Approach,phasing_duration,u_s,u_e,u_b,ts)
%Alex Philpott

%% Phasing dv Calculation
T2 = phasing_duration;
a2 = (T2.*sqrt(u_s)/2/pi)^(2/3);
v = sqrt(u_s.*(2./Earth.a - 1/a2))
dv = sqrt(u_s*(2/Earth.a - 1/a2)) - sqrt(u_s/Earth.a);
dv_phase = dv;

%% Earth departure velocity calculation
v_infin = dv_phase;
vp = sqrt(2*u_e/Launch.P + v_infin^2);
dv1 = vp - Launch.vp;

fprintf('Phasing orbit delta-v: %.3f km/s.\n',dv1);

%% Inclination shift calculation
Earth_pos = [1.501395410244307E+08,-2.768648184673080E+06,-9.023020067359321E+02]; %Earth position on Sept. 22 2017
Earth_r = norm(Earth_pos);

%Hohmann Transfer burn calculation
ra = Earth_r;
rp = Bennu.P;
at = (rp + ra)./2;

vde = sqrt(u_s.*(2./ra - 1./Earth.a));
vdt = sqrt(u_s.*(2./ra - 1./at));
dvd = abs(sqrt(vdt.^2 + vde.^2 - 2*vdt*vde*cosd(abs(Earth.i-Bennu.i)))); %Transfer insertion and inclination change

%Calculate departure burn
v_infin1 = dvd;
vp2 = sqrt(2.*u_e./Launch.P+ v_infin1.^2);
dv1 = vp2-vp;

fprintf('Depart launch orbit with delta-v %.3f km/s. Burn at perigee.\n',dv1);

%% Orbit Simulation
%Positions at Sept 9 2016 07:00:00 UTC
load('parsed_files.mat');
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

y0 = [Sun_x0,Earth_x0,Jupiter_x0,Bennu_x0,OR_x0,Sun_v0,Earth_v0,Jupiter_v0,Bennu_v0,OR_v0]';

t0 = ts(1);
tf = ts(2)*24*3600;
tspan = linspace(t0,tf,ts(2)-ts(1)+1);

options = odeset('AbsTol',1e-10,'RelTol',1e-10);

[Sun_pos,Earth_pos,Jupiter_pos,Bennu_pos,OR_pos,t,Sun_vel,Earth_vel,Jupiter_vel,Bennu_vel,OR_vel] = run_de(tspan,y0,options);

t0 = tf;
tf = ts(3).*24.*3600;
tspan = linspace(t0,tf,ts(3)-ts(2)+1);
OR_dv = (1+.95./100.)*OR_vel(end,:);
fprintf('Implemented DV %f km/s on Day 111.\n',0.95./100.*norm(OR_vel(end,:)));
y0 = [Sun_pos(end,:) Earth_pos(end,:) Jupiter_pos(end,:) Bennu_pos(end,:) OR_pos(end,:) Sun_vel(end,:) ...
    Earth_vel(end,:) Jupiter_vel(end,:) Bennu_vel(end,:) OR_dv];

[Sun_pos2,Earth_pos2,Jupiter_pos2,Bennu_pos2,OR_pos2,t2,Sun_vel2,Earth_vel2,Jupiter_vel2,Bennu_vel2,OR_vel2] = run_de(tspan,y0,options);

Sun_pos = [Sun_pos; Sun_pos2];
Earth_pos = [Earth_pos; Earth_pos2];
Jupiter_pos = [Jupiter_pos; Jupiter_pos2];
Bennu_pos = [Bennu_pos; Bennu_pos2];
OR_pos = [OR_pos; OR_pos2];
t = [t; t2];
OR_vel = [OR_vel; OR_vel2];

end

function [Sun_pos,Earth_pos,Jupiter_pos,Bennu_pos,OR_pos,t,Sun_vel,Earth_vel,Jupiter_vel,Bennu_vel,OR_vel] = run_de(tspan,y0,options)
[t,y] = ode45(@rates,tspan,y0,options);
Sun_pos = y(:,1:3);
Earth_pos = y(:,4:6);
Jupiter_pos = y(:,7:9);
Bennu_pos = y(:,10:12);
OR_pos = y(:,13:15);

Sun_vel = y(:,16:18);
Earth_vel = y(:,19:21);
Jupiter_vel = y(:,22:24);
Bennu_vel = y(:,25:27);
OR_vel = y(:,28:30);
end