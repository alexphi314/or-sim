function dy = rates(~,y)
%Pre-allocating dy
dy = zeros(30,1);
Sun_mass = 1.989e30; %kg
Bennu_mass = 0.00014e15;
Earth_mass = 5.97219e24; %kg
Jupiter_mass = 1.8981e27; %kg
OR_mass = 2110; %kg

%Assigning constants and calculating the current radius
Sun_pos = y(1:3);
Earth_pos = y(4:6);
Jupiter_pos = y(7:9);
Bennu_pos = y(10:12);
OR_pos = y(13:15);

Sun_veloc = y(16:18);
Earth_veloc = y(19:21);
Jupiter_veloc = y(22:24);
Bennu_veloc = y(25:27);
OR_veloc = y(28:30);

X = [Sun_pos';Earth_pos';Jupiter_pos';Bennu_pos';OR_pos'];
M = [Sun_mass,Earth_mass,Jupiter_mass,Bennu_mass,OR_mass];
a = calc_accel(M,X);

% Assigning output
dy(1:3) = [0 0 0];
dy(4:6) = Earth_veloc;
dy(7:9) = Jupiter_veloc;
dy(10:12) = Bennu_veloc;
dy(13:15) = OR_veloc;

dy(16:18) = [0 0 0];
dy(19:21) = a(2,:);
dy(22:24) = a(3,:);
dy(25:27) = a(4,:);
dy(28:30) = a(5,:);
end