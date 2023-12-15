%% Mushy Zone Plot
tm = 0; a = 3; tr = 4;
syms x
sigmoidFunc = 1 / (1 + exp(-a * (x-tm)));
smoothedHeavside = 1/pi * atan((x-tm)/0.1)+0.5;
linear = piecewise(x<-2+tm, 0, (x>-2+tm & x<2+tm), 1/(tr)*(x-tm)+0.5,x>2+tm, 1);

fplot([sigmoidFunc,smoothedHeavside, linear], [-10,10], "LineWidth",2);
ylim([-0.1,1.1]); grid on;
legend("Sigmoid Function", "Smoothed Heaviside", "Linear", Location="northwest");
xlabel("Temperature"); ylabel("Fraction of liquid")


%% Lattice properties
% Structural
syms r
h = 10.0e-3; E = 70e6;
E1bcc = (2*sqrt(2)*pi*E*(r*1e-3)^2/h^2) * (1 + 12*((r*1e-3)^2/h^2)*(sin(pi/4)^2+1)*tan(pi/4)^2) * (sin(pi/4)^2*cos(pi/4));
E3bcc = (8*pi*E*(r*1e-3)^2/h^2) * (1 + 12*((r*1e-3)^2/h^2) * cos(pi/4)^2) * (sin(pi/4)^3 * tan(pi/4)^2);

figure(1); clf(1);
fplot(r, E1bcc/1e6,[0,1.5], 'LineWidth', 2.0); hold on;
fplot(r, E3bcc/1e6,[0,1.5], 'LineWidth', 2.0);
title("Young's Modulus of bcc Unit Cell")
xlabel("Strut radius [mm]"); ylabel("Young's Modulus [MPa]");
legend("1-Direction", "3-Direction", 'Location','northwest');

% Thermal
syms epsilon
phi = pi/4;
h_cell = 10.0;
k_pcm = 0.21;	
k_cell = 160.0;	
rho_pcm = 910.0;
rho_cell = 2650.0;
cp_pcm = 2000.0;
cp_cell = 910.0;
L_pcm = 190000.0;

theta = atan(tan(phi) / sqrt(2));
epsilon = 1 - 2*tan(theta^2* (r^2/h_cell^2) * (pi * (1+4/sin(theta)) -16/3* r/h_cell * (3.137/sin(pi-2*theta) + 4.923/sin(pi/2-theta))));
V4 = (2 * (16 / (3 * sin(pi - 2*theta))) * r^3) - (12 * (sqrt(8) - sqrt(6)) * r^3);
t4 = (V4 * sin(theta) * cos(theta)^2) ^ (1/3);
Lstr = (h_cell / (2 * sin(theta))) - (t4 * sqrt((2 / cos(theta)^2) + (1 / sin(theta)^2)));
R1 = (2 * cos(theta)^2) / (4.277 * sin(theta)* t4 * k_cell);
R2 = Lstr / (pi * r^2 * k_cell);
R3 = R1;
Rstr = 2 * (R1 + R2 + R3);
Rs = Rstr / 4;
k1 = (epsilon * k_pcm) + (1 / (h_cell * Rs));

Lstr = (h_cell / (2 * sin(theta))) - (t4 * sqrt((1 / cos(theta)^2) + (2 / sin(theta)^2)));
R1 = (2 * sin(theta)^2) / (1.903 * cos(theta)* t4 * k_cell);
R2 = Lstr / (pi * r^2 * k_cell);
R3 = R1;
Rstr = 2 * (R1 + R2 + R3);
Rs = Rstr / 4;
k3 = (epsilon * k_pcm) + ((2 * tan(theta)^2) / (h_cell * Rs));

syms porosity
rho = porosity * rho_cell + (1 - porosity) * rho_pcm;
cp = ((porosity * rho_cell) / rho * cp_cell) + (((1 - porosity) * rho_pcm) / rho * cp_pcm);
L = (1 - porosity) * rho_pcm / rho * L_pcm;

figure(2); clf(2); subplot(1,2,1);
fplot(r, k1, [0,1.5],'LineWidth', 2.0); hold on;
fplot(r, k3, [0,1.5],'LineWidth', 2.0);
title("Effective Thermal Conductivity of bcc Unit Cell")
xlabel("Strut radius [mm]"); ylabel("Thermal Conductivity [W/mK]");
legend("1-Direction", "3-Direction", 'Location','northwest');

subplot(1,2,2);
fplot(porosity, L/1e3, [0, 0.4], 'LineWidth', 2.0);
xlabel("Volume Fraction [-]"); ylabel("Latent Heat [kJ]");
title("Effective Latent Heat for bcc Unit Cell");

%% Material properties
syms x
cond = 1000 * piecewise(x < 0, 2.22, x > 0, 0.54);
cp = piecewise(x < 0, 2090, x > 0, 4186) - 2090;
enthalpy = 1e-3 * (x * cp + piecewise(x < 0,0, x > 0, 334000));

fplot([cond, cp, enthalpy],[-10,10], "LineWidth",2.0);

%% Thermal energy storage
syms T

Ts = 25.0;
sensible_cp = 4182;
sensible_rho = 988;
sensible_Q = sensible_rho * sensible_cp * (T - Ts);

Tm = 56.8;
latent_cp = piecewise(T < Tm, 2000, T > Tm, 2150);
latent_rho = piecewise(T < Tm, 910, T > Tm, 790);
latent_dH = piecewise(T < Tm, 0, T > Tm, 190000);
latent_Q = latent_rho * latent_cp * (T - Ts) + latent_rho * latent_dH;

fplot(1e-6*[sensible_Q, latent_Q], [Ts, 90], "LineWidth",2.0)
xlabel("Temperature [°C]"); ylabel("Stored enerrgy [MJ]")
legend(["Sensible storage", "Latent heat storage"], "Location","northwest")


%% PCM Properties Selection
colours = [
    0 0.4470 0.7410;
    0.8500 0.3250 0.0980;
    0.9290 0.6940 0.1250;
    0.4940 0.1840 0.5560;
    0.4660 0.6740 0.1880;
    0.3010 0.7450 0.9330;
    0.6350 0.0780 0.1840;
    0 0.4470 0.7410;
    0.8500 0.3250 0.0980;
    0.9290 0.6940 0.1250;
    0.4940 0.1840 0.5560;
    0.4660 0.6740 0.1880;
    0.3010 0.7450 0.9330;
    0.6350 0.0780 0.1840;];

materials = [
    "Water-Salt Solutions";
    "Water";
    "Clathrates";
    "Paraffins";
    "Salt Hydrates";
    "Sugar Alcohols";
    "Nitrates";
    "Hydroxides";
    "Chlorides";
    "Carbonates";
    "Flourides"];

properties = [
    -100, 0, 200, 300;
    -2, 2, 328, 332;
    -50, 0, 200, 300;
    -20, 100, 150, 250;
    -20, 80, 200, 600;
    20, 450, 200, 450;
    120, 300, 200, 700;
    150, 400, 500, 700;
    350, 750, 550, 800;
    400, 800, 600, 1000;
    700, 900, 950, 1100];

figure(1); clf(1);
for i = 1:length(materials)
    [major, minor, center] = getMajorMinorCenter(properties(i,1),properties(i,2),properties(i,3),properties(i,4));
    plotEllipse(major, minor, center, colours(i,:));
end
for i = 1:length(materials)
    [major, minor, center] = getMajorMinorCenter(properties(i,1),properties(i,2),properties(i,3),properties(i,4));
    t = text(center(1), center(2), materials(i));
    t.Position(1) = center(1) - 0.5 * t.Extent(3);
end
xlabel("Temperature [°C]"); ylabel("Melting Latent Heat");

function [major, minor, center] = getMajorMinorCenter(minX, maxX, minY, maxY)
    major = 0.5 * (maxX - minX);
    minor = 0.5 * (maxY - minY);
    center = 0.5 * [maxX + minX, maxY + minY];
end

function plotEllipse(major, minor, center, c)
    alpha = linspace(0,2*pi,100)';
    p = [major * cos(alpha), minor * sin(alpha)];
    patch(p(:,1) + center(1), p(:,2) + center(2), c, 'FaceAlpha',0.5); hold on;
end



