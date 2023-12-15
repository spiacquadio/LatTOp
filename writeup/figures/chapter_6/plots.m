%% Case 1 =================================================================
%% MF Compliance
c = readmatrix("complianceCase2.csv");
x = linspace(0, 1, length(c));
structures = (x == 0 | x == 0.25 | x == 0.5 | x == 0.75 | x == 1 );

strings = { "Weight = 0.0", "Weight = 0.25", "Weight = 0.5", "Weight = 0.75", 'Weight = 1.0' };

f = figure(1); clf(1);
plot(c(:,1), c(:,2), 'LineWidth',2.0); hold on;
plot(c(structures,1), c(structures,2),'b*');
text(c(structures,1), c(structures,2), strings);
plot(utopiaPoint(1),utopiaPoint(2),'m*');
title("Pareto Front");
xlabel("Thermal Compliance");
ylabel("Structural Compliance");

utopiaPoint = [312.6,0.01]; wopt = 0.49194024058363417;
copt = [312.6229473221757, 0.011612970464663936];
c1range = max(c(:,1)) - min(c(:,1));
c2range = max(c(:,2)) - min(c(:,2));
gopt = goalFunction(copt, utopiaPoint, c1range,c2range);
g = goalFunction(c, utopiaPoint, c1range,c2range);
figure(2); clf(2);
plot(x, g, "LineWidth", 2.0); hold on;
plot(wopt, gopt, "r*");
xlabel("Weight");
ylabel("L_3 Error norm")


%% Case 2 ================================================================
%% MF Compliance
c = readmatrix("complianceCase1.csv");
x = linspace(0, 1, length(c));
structures = (x == 0 | x == 0.25 | x == 0.5 | x == 0.75 | x == 1 );

strings = { "Weight = 0.0", "Weight = 0.25", "Weight = 0.5", "Weight = 0.75", 'Weight = 1.0' };

f = figure(1); clf(1);
plot(c(:,1), c(:,2), 'LineWidth',2.0); hold on;
plot(c(structures,1), c(structures,2),'b*');
text(c(structures,1), c(structures,2), strings);
plot(utopiaPoint(1),utopiaPoint(2),'m*');
title("Pareto Front");
xlabel("Thermal Compliance");
ylabel("Structural Compliance");
% xlim([312.59,312.85])

utopiaPoint = [1012.66,0.08]; wopt = 0.36446256247444675;
copt = [1012.6608717696596, 0.08145351603379186];
c1range = max(c(:,1)) - min(c(:,1));
c2range = max(c(:,2)) - min(c(:,2));
gopt = goalFunction(copt, utopiaPoint, c1range,c2range);
g = goalFunction(c, utopiaPoint, c1range,c2range);
figure(2); clf(2);
plot(x, g, "LineWidth", 2.0); hold on;
plot(wopt, gopt, "r*");
xlabel("Weight");
ylabel("L_3 Error norm")


%% Case 3 ================================================================
%% MF Compliance
c = readmatrix("complianceCase3.csv");
x = linspace(0, 1, length(c));
structures = (x == 0 | x == 0.25 | x == 0.5 | x == 0.75 | x == 1 );

strings = { "Weight = 0.0", "Weight = 0.25", "Weight = 0.5", "Weight = 0.75", 'Weight = 1.0' };

utopiaPoint = [312.6,0.08]; wopt = 0.5818856709054164;
copt = [312.6081771824148, 0.08118559675891564];
c1range = max(c(:,1)) - min(c(:,1));
c2range = max(c(:,2)) - min(c(:,2));
gopt = goalFunction(copt, utopiaPoint, c1range,c2range);
g = goalFunction(c, utopiaPoint, c1range,c2range);

f = figure(1); clf(1);
plot(c(:,1), c(:,2), 'LineWidth',2.0); hold on;
plot(c(structures,1), c(structures,2),'b*');
text(c(structures,1), c(structures,2), strings);
plot(utopiaPoint(1),utopiaPoint(2),'m*');
title("Pareto Front");
xlabel("Thermal Compliance");
ylabel("Structural Compliance");

figure(2); clf(2);
plot(x, g, "LineWidth", 2.0); hold on;
plot(wopt, gopt, "r*");
xlabel("Weight");
ylabel("L_3 Error norm")


%% Projection Idea
syms x
k = 0.1;

a = -0.5 + 0.5*sqrt(1+4*k);
b = k / (1+a);
f = k/(x+a) - b;

u = (utopiaPoint - [min(c(:,1)), min(c(:,2))]) ./ [c1range,c2range];
g = matlabFunction(norm([f-u(2),x-u(1)],3));
ff = matlabFunction(f);
gOpt = fminsearch(g,0.5);

ds = sqrt(1 + gradient(f,x)^2);
sFull = eval(int(ds,0,1));
w = 1 - (eval(int(ds,0,gOpt)) / sFull);

gGuess = getInitialGuess(c(1,:)./[c1range,c2range],c(end,:)./[c1range,c2range],utopiaPoint./[c1range,c2range]);

figure(3); clf(3);
fplot(f,[0,1],"LineWidth",2.0); hold on;
% plot([0,1],[1,0],"Color",[0, 0.4470, 0.7410],"LineStyle","--","LineWidth",2.0);
plot(u(1),u(2),"r*")
plot(gOpt,ff(gOpt), "b*");
title("Utopia Point Projected onto Reciprocal Function");
ylim([0,1]);

%% Helper Functions
function err = goalFunction(c, utopia, c1range, c2range)
    err = vecnorm([(c(:,1) - utopia(1))/c1range, (c(:,2) - utopia(2))/c2range], 3, 2);
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







