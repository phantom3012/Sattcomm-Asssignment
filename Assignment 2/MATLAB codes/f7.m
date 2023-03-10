clear
close all
clc

%%%%%%%%%%%%%%%%%%%%
%Required constants%
%%%%%%%%%%%%%%%%%%%%

a = 6589116;
e = 0.007589;
i = 32.54;
Omega = 235.2;
omega = 181.2;
M0 = 228.5;
T0 = 2437716.11642;

GMe = 3.986004415 *10^(14);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Taking required input from the user%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp("Enter the precise date and time needed to find latitiude and longitude of the Friendship 7 spacecraft")

date = input("Enter the required date: ");
month = input ("Enter the required month(in number format): ");
Y = input ("Enter the required year: ");
hh = input("Enter the required hour of the day (24h format): ");
mm = input("Enter the required minute of the hour: ");
ss = input("Enter the required second of the minute: ");

d = ((hh/24)+(mm/1440)+(ss/86400));
D = date + d;
if month<=2
    month = month+12;
    Y = Y -1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of the Julian day
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = floor(Y/100);
B = 2 - A + floor(A/4);

t = (floor(365.25*(Y+4716))+floor(30.6001*(month+1))+D+B-1524.5);
% t = 2459580;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating required orbital parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = (86400/(2*pi))*sqrt(GMe/a^3);
M = deg2rad(M0) + 2*pi*n*(t-T0);

syms E_sym
eqn = M == E_sym - e*sin(E_sym);
E_r = solve(eqn,E_sym);
E = rad2deg(double(E_r));
f = double(2*atand(sqrt((1+e)/(1-e)) * tan (E/2)));
u = omega + f;

r = a*(1-e*cosd(E));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Converting orbital parameters to Latitude and Longitude
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%cartesian conversion

x = r*((cosd(u)*cosd(Omega))-(sind(u)*sind(Omega)*cosd(i)));
y = r*((cosd(u)*sind(Omega))+(sind(u)*cosd(Omega)*cosd(i)));
z = r*sind(u)*sind(i);

Alpha = atan2d(y,x);
Delta = asind(z/r);

T = ((t-2451545)/36525);

GST_d = (280.46061837 + (360.98564736629*(t-2451545)) + (0.000387933*T^2) - ((T^3)/38710000));
GST = mod(GST_d,360); 

L = wrapTo180(alpha - GST);
psi = Delta;

psi_s = angl2str(psi,"ns","degrees2dms");
L_s = angl2str(L,"ew","degrees2dms");

sprintf("The latitude of the spacecraft at this time is %s ",psi_s)
sprintf("The longitude of the spacecraft at this time is: %s\n",L_s)