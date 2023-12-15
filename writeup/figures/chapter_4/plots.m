%% Phase Change Simulations
dt = 30.0;
pcGeneric = readmatrix("PhaseChangeDataGeneric.csv");
pcLattice0to20 = readmatrix("PhaseChangeDataLattice0to20.csv");
pcLattice5to20 = readmatrix("PhaseChangeDataLattice5to20.csv");
pcLattice10to20 = readmatrix("PhaseChangeDataLattice10to20.csv");

plot(pcLattice5to20(:,2),[pcLattice5to20(:,1),pcLattice0to20(:,1),pcLattice10to20(:,1),pcGeneric(:,1)], ...
    "LineWidth",2.0);
xlim([0,2800]); ylim([295,400]);
title("Phase Change Comparison");
legend(["5% to 20%", "0% to 20%", "10% to 20%", "Generic Properites"], "Location","northwest");
ylabel("Temperature [K]"); xlabel("Time [s]")

time = 3500;
idx = find(abs(pcGeneric(:,2) - 3500) == min(abs(pcGeneric(:,2)-3500)));
relError = 100 .* [(pcLattice0to20(idx,1) - pcLattice5to20(idx,1))/pcLattice5to20(idx,1) ...
    (pcGeneric(idx,1) - pcLattice5to20(idx,1))/pcLattice5to20(idx,1)]

%%
Tm = 273.15 + 56.8;
lfGeneric = 1 ./ (1 + exp(-0.5*(pcGeneric(:,1) - Tm)));
lf0to20 = 1 ./ (1 + exp(-0.5*(pcLattice0to20(:,1) - Tm)));
lf5to20 = 1 ./ (1 + exp(-0.5*(pcLattice5to20(:,1) - Tm)));
lf10to20 = 1 ./ (1 + exp(-0.5*(pcLattice10to20(:,1) - Tm)));

plot(pcLattice5to20(:,2),[lf5to20,lf0to20,lf10to20,lfGeneric],"LineWidth",2.0);
xlim([0,600])
title("Phase Change Comparison");
legend(["5% to 20%", "0% to 20%", "10% to 20%", "Generic Properites"], "Location","northwest");
ylabel("Temperature [K]"); xlabel("Time [s]")


%% Phase Change Optimums PCsims
dt = 30.0;
pcCondDom = readmatrix("PhaseChangeConductionDominant.csv");
pcLatentDom = readmatrix("PhaseChangeLatentHeatDominant.csv");
pcIncom = readmatrix("PhaseChangeIncompleteMelting.csv");
pcPureCond = readmatrix("PhaseChangePureConduction.csv");

plot(pcCondDom(:,2),[pcCondDom(:,1),pcLatentDom(:,1),pcIncom(:,1),pcPureCond(:,1)], ...
    "LineWidth",2.0);
xlim([0,3000]); ylim([295,400]);
title("Phase Change Comparison");
legend(["Conduction Dominant", "Latent Heat Dominant", "Incomplete Melting", "Pure Conduction"], "Location","northwest");
ylabel("Temperature [K]"); xlabel("Time [s]")