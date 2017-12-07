function [theta_b] = calc_ta_bennu(Bennu,dt,theta_b)
%Alex Philpott
f_prime = @(x) 1 - Bennu.e.*cos(x);
Mb = Bennu.n.*dt;
fhb = @(x) x - Bennu.e.*sin(x) - Mb;
Eb = newton_method(fhb,f_prime,theta_b,1e-6);
theta_b_new = 2*atan(sqrt((1+Bennu.e)/(1-Bennu.e))*tan(Eb/2));
theta_b = theta_b_new + theta_b;
end