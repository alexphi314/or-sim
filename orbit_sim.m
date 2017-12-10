function [Sun_pos,Earth_pos,Jupiter_pos,Bennu_pos,OR_pos,t,OR_vel] = orbit_sim(Bennu,Earth,Launch,Survey,Approach,phasing_duration,u_s,u_e,u_b,ts)
%Alex Philpott

%% Phasing dv Calculation
T2 = phasing_duration;
a2 = (T2.*sqrt(u_s)/2/pi)^(2/3);
v = sqrt(u_s.*(2./Earth.a - 1/a2));
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

%% Simulating
options = odeset('AbsTol',1e-10,'RelTol',1e-10);
y0 = [Sun_x0,Earth_x0,Jupiter_x0,Bennu_x0,OR_x0,Sun_v0,Earth_v0,Jupiter_v0,Bennu_v0,OR_v0]';
Sun_pos = [];
Sun_pos(1,:) = Sun_x0;
Earth_pos = [];
Earth_pos(1,:) = Earth_x0;
Jupiter_pos = [];
Jupiter_pos(1,:) = Jupiter_x0;
Bennu_pos = [];
Bennu_pos(1,:) = Bennu_x0;
OR_pos = [];
OR_pos(1,:) = OR_x0;
Sun_vel = [];
Sun_vel(1,:) = Sun_v0;
Earth_vel = [];
Earth_vel(1,:) = Earth_v0;
Jupiter_vel = [];
Jupiter_vel(1,:) = Jupiter_v0;
Bennu_vel = [];
Bennu_vel(1,:) = Bennu_v0;
OR_vel = [];
OR_vel(1,:) = OR_v0;
t = [];
t(1) = 0;
for k = 2:7
    t0 = ts(k-1)*24*3600;
    tf = ts(k)*24*3600;
    tspan = linspace(t0,tf,ts(k)-ts(k-1)+1);

    [Sun_pos2,Earth_pos2,Jupiter_pos2,Bennu_pos2,OR_pos2,t2,Sun_vel2,Earth_vel2,Jupiter_vel2,Bennu_vel2,OR_vel2]...
        = run_de(tspan,y0,options);
    
    Sun_pos = [Sun_pos; Sun_pos2(2:end,:)];
    Sun_vel = [Sun_vel; Sun_vel2(2:end,:)];
    Earth_pos = [Earth_pos; Earth_pos2(2:end,:)];
    Earth_vel = [Earth_vel; Earth_vel2(2:end,:)];
    Jupiter_pos = [Jupiter_pos; Jupiter_pos2(2:end,:)];
    Jupiter_vel = [Jupiter_vel; Jupiter_vel2(2:end,:)];
    Bennu_pos = [Bennu_pos; Bennu_pos2(2:end,:)];
    Bennu_vel = [Bennu_vel; Bennu_vel2(2:end,:)];
    OR_pos = [OR_pos; OR_pos2(2:end,:)];
    OR_vel = [OR_vel; OR_vel2(2:end,:)];
    t = [t; t2(2:end,:)];
    
    y0 = [Sun_pos(end,:) Earth_pos(end,:) Jupiter_pos(end,:) Bennu_pos(end,:) OR_pos(end,:) Sun_vel(end,:) ...
    Earth_vel(end,:) Jupiter_vel(end,:) Bennu_vel(end,:) OR_vel(end,:)];

    if ts(k) == 111
        OR_dv = (1+.95./100.).*OR_vel(end,:);
        y0(28:30) = OR_dv;
        fprintf('Implemented DV %f km/s on Day %i.\n',0.95./100.*norm(OR_vel(end,:)),ts(k));
    end
    
    if ts(k) == 378
        y0(13:15) = E_OR(ts(3),:); %Moving OR position
        y0(28:30) = VE_OR(ts(3),:); %Moving OR velocity
        y0(4:6) = E_Earth(ts(3),:); %Moving Earth position
        y0(19:21) = VE_Earth(ts(3),:); %Moving Earth velocity
    end
    
    if ts(k) == 754
        %OR_dv = (1+1.394/100).*OR_vel(end,:);
        y0(28:30) = VE_OR(754,:);
        y0(13:15) = E_OR(754,:);
        y0(10:12) = E_Bennu(754,:);
        y0(25:27) = VE_Bennu(754,:);
        %fprintf('Implemented DV %.3f km/s on Day %i.\n',1.34/100*norm(OR_vel(end,:)),ts(k));
    end
    
    if ts(k) == 768
        OR_dv = (1+0.739/100).*OR_vel(end,:)
        y0(28:30) = OR_dv;
        VE_OR(768,:)
        fprintf('Implemented DV %.3f km/s on Day %i.\n',0.739/100*norm(OR_vel(end,:)),ts(k));
    end
end
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