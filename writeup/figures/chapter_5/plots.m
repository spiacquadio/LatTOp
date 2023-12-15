%% Multifunctional Phase-Change Simulations
pc1to0 = readmatrix("PC_MF_1to0.csv");
pc3to1 = readmatrix("PC_MF_3to1.csv");
pc1to1 = readmatrix("PC_MF_1to1.csv");
pc1to3 = readmatrix("PC_MF_1to3.csv");
pc0to1 = readmatrix("PC_MF_0to1.csv");

plot(pc1to0(:,2),[pc1to0(:,1),pc3to1(:,1),pc1to1(:,1),pc1to3(:,1),pc0to1(:,1)], "LineWidth",2.0);
xlim([0,3000]); ylim([275,425]);
title("Multi-Functional Structure Phase-Change Results")
legend(["Thermal to Solid ratio 1:0", "Thermal to Solid ratio 3:1", ...
    "Thermal to Solid ratio 1:1", "Thermal to Solid ratio 1:3", ...
    "Thermal to Solid ratio 0:1"], "Location","northwest");

%% MF Compliance
c = readmatrix("compliance.csv");
x = linspace(0, 1, length(c));
structures = (x == 0 | x == 0.25 | x == 0.5 | x == 0.75 | x == 1 );

strings = { "Weight = 0.0", "Weight = 0.25", "Weight = 0.5", "Weight = 0.75", 'Weight = 1.0' };

f = figure(1); clf(1);
plot(c(:,1), c(:,2), 'LineWidth',2.0); hold on;
plot(c(structures,1), c(structures,2),'b*');
text(c(structures,1), c(structures,2), strings);
title("Pareto Front");
xlabel("Thermal Compliance");
ylabel("Structural Compliance");


%% MF Optimum
utopia = [1012.66, 0.011];
cOpt = [1012.6643974234428, 0.01170652390849104]; wOpt = 0.63; gOpt = 0.07123166053283343;
% utopia = [1012.64, 0.014];
% cOpt = [1012.6478748412313, 0.014556708347830391]; wOpt = 0.85399749022446; gOpt = 0.09073288366725593;
c = readmatrix("compliance.csv");

x = linspace(0, 1, length(c));
c1range = max(c(:,1)) - min(c(:,1));
c2range = max(c(:,2)) - min(c(:,2));
goal = goalFunction(c, utopia, c1range, c2range);

figure(1); clf(1); subplot(1,2,1);
plot(c(:,1), c(:,2), 'LineWidth',2.0); hold on;
plot(utopia(1),utopia(2), 'm*');
plot(cOpt(1),cOpt(2), 'r*');
title("Pareto Front");
xlabel("Thermal Compliance");
ylabel("Structural Compliance");

subplot(1,2,2);
plot(x, goal, 'LineWidth', 2.0); hold on;
plot(wOpt,gOpt, 'r*');
title("Error Measure");
xlabel("Weight"); ylabel('L_3 Error');


%% Constraint Method Limitation
n = 2000;
x = linspace(0, 1, n)';
utopia = [-600, 20];

c = zeros(n,2);
for i = 1:n
    c(i,:) = [f1(x(i)), f2(x(i))];
end

c2range = c(end,2) - c(1,2);
c1range = c(1,1) - c(end,1);
goal = goalFunction(c, utopia, c1range, c2range);

c1c = abs(c(:,1) - utopia(1));
c2c = abs(c(:,2) - utopia(2));
cci = (c1c == min(c1c)) | (c2c == min(c2c));

figure(1); clf(1); subplot(1,2,1);
plot(c(:,1), c(:,2), 'LineWidth',2.0); hold on;
plot(utopia(1), utopia(2), "m*");
plot(c(goal==min(goal),1), c(goal==min(goal),2), 'r*')
plot(c(cci,1),c(cci,2),'b*');
title("Pareto Front");
xlim([-800,100]);

subplot(1,2,2);
plot(x, goal, 'LineWidth', 2.0); hold on;
title("L_2 Error");
plot(x(goal==min(goal)), min(goal), 'r*');
plot(x(cci), goal(cci), 'b*');


%% Quadratic Fit Limitations
n = 2000;
x = linspace(0, 1, n)';
utopia = [-900, 34];

c = zeros(n,2);
for i = 1:n
    c(i,:) = [f1(x(i)), f2(x(i))];
end

c2range = c(end,2) - c(1,2);
c1range = c(1,1) - c(end,1);
goal = goalFunction(c, utopia, c1range, c2range);

figure(1); clf(1); subplot(1,2,1);
plot(c(:,1), c(:,2), 'LineWidth',2.0); hold on;
plot(utopia(1), utopia(2), "m*");
plot(c(goal==min(goal),1), c(goal==min(goal),2), 'r*')
title("Pareto Front");
% xlim([-800,100]);

subplot(1,2,2);
plot(x, goal, 'LineWidth', 2.0); hold on;
title("L_2 Error");
plot(x(goal==min(goal)), min(goal), 'r*');


%% MF Approach
utopia = [1012.66, 0.011];
cOpt = [1012.6643974234428, 0.01170652390849104]; wOpt = 0.63; gOpt = 0.07123166053283343;
% utopia = [1012.64, 0.014];
% cOpt = [1012.6478748412313, 0.014556708347830391]; wOpt = 0.85399749022446; gOpt = 0.09073288366725593;
c = readmatrix("compliance.csv");

x = linspace(0, 1, length(c));
c1range = max(c(:,1)) - min(c(:,1));
c2range = max(c(:,2)) - min(c(:,2));
goal = goalFunction(c, utopia, c1range, c2range);
cTmax = max(c(:,1)); cTmin = min(c(:,1));
cSmax = max(c(:,2)); cSmin = min(c(:,2));
[w,~,~] = getInitialGuess([cTmax,cSmin]./[c1range, c2range], ...
    [cTmin,cSmax]./[c1range, c2range], utopia./[c1range, c2range]);

idxClosest = find(abs(w-x) == min(abs(w-x)));
v = [cTmax,cSmin] - [cTmin,cSmax];
pg = -v*w + [cTmax,cSmin];

figure(1); clf(1); subplot(1,2,1);
plot(c(:,1), c(:,2), "--", 'LineWidth',1.0); hold on;
plot(c([1,end],1), c([1,end],2), "LineWidth",2.0,"Color",[0 0.4470 0.7410]);
plot(utopia(1),utopia(2), 'm*');
plot(cOpt(1),cOpt(2), 'r*');
plot(pg(1),pg(2),"b*")
title("Pareto Front"); xlabel("Thermal Compliance"); ylabel("Structural Compliance");
legend("True Pareto Front","Construction Line Segment", Location="northeast")

subplot(1,2,2);
plot(x, goal, 'LineWidth', 2.0); hold on;
plot(wOpt,gOpt, 'r*');
plot(x(idxClosest), goal(idxClosest),"b*")
title("Error Measure");
xlabel("Weight"); ylabel('L_3 Error');

%% Convergence
utopia = [1012.66, 0.011];
cOpt = [1012.6643974234428, 0.01170652390849104]; wOpt = 0.63; gOpt = 0.07123166053283343;
% utopia = [1012.64, 0.014];
% cOpt = [1012.6478748412313, 0.014556708347830391]; wOpt = 0.85399749022446; gOpt = 0.09073288366725593;
c = readmatrix("compliance.csv");

x = linspace(0, 1, length(c));
c1range = max(c(:,1)) - min(c(:,1));
c2range = max(c(:,2)) - min(c(:,2));
goal = goalFunction(c, utopia, c1range, c2range);

data = [0.0       0.422914  0.471136  0.573958   0.60441    0.60441    0.60441    0.60441    0.60441    0.628878   0.628878;
 0.422914  0.471136  0.573958  0.60441    0.63453    0.634017   0.630797   0.630117   0.62908    0.62908    0.62908;
 1.0       1.0       1.0       1.0        1.0        0.63453    0.634017   0.630797   0.630117   0.630117   0.629263;
 0.846026  0.25481   0.197216  0.099629   0.0804567  0.0804567  0.0804567  0.0804567  0.0804567  0.0741883  0.0741883;
 0.25481   0.197216  0.099629  0.0804567  0.0743916  0.0743599  0.0742207  0.0742014  0.0741857  0.0741857  0.0741857;
 1.00141   1.00141   1.00141   1.00141    1.00141    0.0743916  0.0743599  0.0742207  0.0742014  0.0742014  0.0741925]';

pa = data(:,1); pb = data(:,2); pc = data(:,3);
ga = data(:,4); gb = data(:,5); gc = data(:,6);

figure(1); clf(1);
for i = 1:3
    subplot(2,3,i);
    plot(x, goal, "LineWidth", 2.0); hold on;
    plotParabola([pa(i),ga(i)], [pb(i),gb(i)], [pc(i),gc(i)]);
    ylim([0,1.1]); xlim([0,1]);
    xlabel("Weight"); ylabel('L_3 Error');
    title(sprintf("Iteration: %i", i));
end
for i = 4:6
    subplot(2,3,i);
    plot(x, goal, "LineWidth", 2.0); hold on;
    plotParabola([pa(i),ga(i)], [pb(i),gb(i)], [pc(i),gc(i)]);
    ylim([0,0.25]); xlim([0.5,0.75]);
    xlabel("Weight"); ylabel('L_3 Error');
    title(sprintf("Iteration: %i", i));
end


%% MF Compliance
c = readmatrix("compliance_pc_added.csv");
xtemp = linspace(0, 1, length(c)-2);
x = zeros(13,1);
x(1:3) = xtemp(1:3); x(4) = 0.25; x(5:9) = xtemp(4:8); x(10) = 0.75; x(11:end) = xtemp(9:11);
structures = (x == 0 | x == 0.25 | x == 0.5 | x == 0.75 | x == 1 );

strings = { "Weight = 0.0", "Weight = 0.25", "Weight = 0.5", "Weight = 0.75", 'Weight = 1.0' };

f = figure(1); clf(1);
plot(c(:,1), c(:,2), 'LineWidth',2.0); hold on;
plot(c(structures,1), c(structures,2),'b*');
text(c(structures,1), c(structures,2), strings);
title("Pareto Front");
xlabel("Phase Change Compliance");
ylabel("Structural Compliance");

utopia = [212500,0.012];
cOpt = [213022.65866937238, 0.011989116404131125]; wOpt = 0.5604513990050843; gOpt = 0.087907570508819;

c1range = max(c(:,1)) - min(c(:,1));
c2range = max(c(:,2)) - min(c(:,2));
goal = goalFunction(c, utopia, c1range, c2range);

figure(2); clf(2); subplot(1,2,1);
plot(c(:,1), c(:,2), 'LineWidth',2.0); hold on;
plot(utopia(1),utopia(2), 'm*');
plot(cOpt(1),cOpt(2), 'r*');
title("Pareto Front");
xlabel("Thermal Compliance");
ylabel("Structural Compliance");

subplot(1,2,2);
plot(x, goal, 'LineWidth', 2.0); hold on;
plot(wOpt,gOpt, 'r*');
title("Error Measure");
xlabel("Weight"); ylabel('L_3 Error');


%% Phase Change Temperature
pc0 = readmatrix("PC_MF_w00_coarsest.csv");
pc025 = readmatrix("PC_MF_w025_coarsest.csv");
pc05 = readmatrix("PC_MF_w05_coarsest.csv");
pc075 = readmatrix("PC_MF_w075_coarsest.csv");
pc1 = readmatrix("PC_MF_w1_coarsest.csv");
pcpc0 = readmatrix("PC_MFPC_w00_coarsest.csv");
pcpc025 = readmatrix("PC_MFPC_w025_coarsest.csv");
pcpc05 = readmatrix("PC_MFPC_w05_coarsest.csv");
pcpc075 = readmatrix("PC_MFPC_w075_coarsest.csv");
pcpc1 = readmatrix("PC_MFPC_w1_coarsest.csv");

figure(1); clf(1);
plot(pc0(:,2),[pc0(:,1),pc025(:,1),pc05(:,1),pc075(:,1),pc1(:,1)], "LineWidth",2.0); hold on;
xlim([0,3000]); ylim([275,425]);
title("Multi-Functional Structure Phase-Change Results")

set(gca,'ColorOrderIndex',1)
plot(pc0(:,2),[pcpc0(:,1),pcpc025(:,1),pcpc05(:,1),pcpc075(:,1),pcpc1(:,1)], "LineWidth",2.0,"LineStyle","--");
legend(["TS 1:0", "TS 3:1", ...
    "TS 1:1", "TS 1:3", ...
    "TS 0:1","PS 1:0", ...
    "PS 3:1", "PS 1:1", ...
    "PS 1:3", "PS 0:1"], "Location","northwest");


%% Helper Functions
function err = goalFunction(c, utopia, c1range, c2range)
    err = vecnorm([(c(:,1) - utopia(1))/c1range, (c(:,2) - utopia(2))/c2range], 3, 2);
end

function c = f1(x)
    cmin = 100;
    cmax = 1000;
    c = (cmax - cmin) .* -x.^.3 + cmin;
end

function c = f2(x)
    cmin = 10;
    cmax = 1;
    c = (cmax - cmin) * -exp(x) + cmin;
end

function [p, vec, veclength] = getInitialGuess(p1, p2, px)
	vec = p1 - p2;
	veclength = sqrt(vec(1)^2 + vec(2)^2);
	vec = vec ./ veclength;

	vecpx = p1 - px;
	distance = vecpx(1)*vec(1) + vecpx(2)*vec(2);

	p = distance / veclength;
end

function [a,b,c] = fitParabola(p1,p2,p3)
    A = [p1(1)^2 p1(1) 1;
         p2(1)^2 p2(1) 1;
         p3(1)^2 p3(1) 1];
    b = [p1(2); p2(2); p3(2)];

    x = A \ b;
    a = x(1); b = x(2); c = x(3);
end

function plotParabola(p1,p2,p3)
    [a,b,c] = fitParabola(p1,p2,p3);
    p = [p1;p2;p3];
    % xrange = [min(p(:,1)), max(p(:,1))];
    xrange = [0,1.001];

    syms x
    f = a*x^2 + b*x + c;
    
    fplot(f, xrange,"r--"); hold on;
    plot(p(:,1),p(:,2),"b*");
end







