clear; 

%% Apparent Heat Capacity Method
syms T;
Tm = 56.8;
Ts = 25.0;
L = 190000;
latent_cp = piecewise(T < Tm, 2000, T > Tm, 2150);
latent_rho = piecewise(T < Tm, 910, T > Tm, 790);
latent_dH = piecewise(T < Tm, 0, T > Tm, 190000);
latent_Q = latent_rho * latent_cp * (T - Ts) + latent_rho * latent_dH;

a = 1.25;
sigmoidFunc = 1 / (1 + exp(-a * (T-Tm)));
dsdx = diff(sigmoidFunc, T);

capp = latent_cp + L*dsdx;
fplot(capp, [25, 100], 'LineWidth',2.0);
title("Apparent Heat Capacity");
xlabel("Temperature [°C]");
ylabel("c_{app} [J/kgK]");