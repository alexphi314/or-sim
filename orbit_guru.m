function [a,e,i,RAAN,wp,h,theta] = orbit_guru(r,v,mew)
%Given r and v (in ECI coordinate frame) and mew, return the Keplerian orbit
%elements. Returned in order: a;e;i;RAAN;wp;h;theta

%Quotes (to give you orbit advice when you run the orbit_guru function)
quotes = {'Let go your earthly tether. Enter the void. Empty and become an orbiting object';
 'New orbits cannot exist without first the destruction of the old orbit';
 'Instinct is a lie, told by a fearful body, hoping to be in orbit';
 'At times it is better to know when to use a parabolic escape orbit than a hyperbolic escape orbit';
 'Shoot for the moon! If you miss you will be in an hyperbolic escape orbit';
 'A satellite in LEO is worth three in GEO!';
 'Others enjoy your orbit';
 'If you never expect to achieve orbit you will never be disappointed';
 'A comet of great importance may reach you any day now';
 'Let your orbit wander';
 'One is not in a hyperbolic orbit, does not mean one is in a parabolic orbit';
 'If your orbit is still in one piece, buy lotto';
 'The cure for grief is circular orbits';
 'Judge each day not by the manuevers you plan but the collisions you avoid';
 'If you never give up on parabolic orbits, they will never give up on you';
 'All escape orbits happen for a reason';
 'A new satellite is on the horizon';
 'Maybe its an orbit, maybe its maybelline';};
%Picking a random quote from the above cell
l = size(quotes);
num = randi(l(1),1);
quote = quotes{num,1};

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

disp(quote);
end


