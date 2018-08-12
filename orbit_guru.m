function [a,e,i,RAAN,wp,h,theta] = orbit_guru(r,v,mew)
%Given r and v (in ECI coordinate frame) and mew, return the Keplerian orbit
%elements. Returned in order: a;e;i;RAAN;wp;h;theta

%Magnitudes of vectors
magn_r = norm(r);
magn_v = norm(v);

%Find semi-major axis
E = magn_v.^2./2 - mew./magn_r;
a = -mew./2./E;

%Find angular momentum
h = cross(r,v);
magn_h = norm(h);

%Find eccentricity
e = cross(v,h)./mew - r./magn_r;
magn_e = norm(e);

%Find inclination
i = acos(dot(h,[0 0 1])./magn_h);

%Define n
n = cross([0 0 1],h);
magn_n = norm(n);

%Find RAAN
RAAN = acos(dot(n,[1 0 0])./magn_n);
if n(2) < 0
    RAAN = 2.*pi - RAAN;
end

%Find wp
wp = acos(dot(n,e)./magn_n./magn_e);
if e(3) < 0
    wp = 2.*pi - wp;
end

%Find theta
theta = acos(dot(r,e)./magn_e./magn_r);
if dot(r,v) < 0
    theta = 2.*pi - theta;
end

end


