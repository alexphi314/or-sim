function [tolerances,theta_e,theta_b] = sim_angle(Earth,Bennu,f_prime,theta_e,theta_b,dt,phi,time)
%Alex Philpott
f_primeb = f_prime;
E_prevb = pi/4;
tolerances = [];
theta_eprev = theta_e;
theta_bprev = theta_b;
t = 0;
if phi == 0
    var = time;
else
    var = 10*dt*365.25;
end
while t < var
    Mb = Bennu.n.*dt;
    
    fhb = @(x) x - Bennu.e.*sin(x) - Mb;
    
    Eb = newton_method(fhb,f_primeb,E_prevb,1e-6);
    
    theta_b = 2*atan(sqrt((1+Bennu.e)/(1-Bennu.e))*tan(Eb/2)) + theta_bprev;
    theta_e = Earth.n.*dt + theta_eprev;
    
    if theta_e > 2*pi && theta_eprev < 2*pi
        theta_e = theta_e - 2*pi;
    end
    
    if theta_b > 2*pi && theta_bprev < 2*pi
        theta_b = theta_b - 2*pi;
    end
    
    diff = theta_e - theta_b;
    tolerances = [tolerances abs(diff - phi)/abs(phi)];
    
    E_prevb = Eb;
    
    theta_bprev = theta_b;
    theta_eprev = theta_e;
    
    t = t + dt;
end