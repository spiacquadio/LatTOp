%% Read Data from file
data = readmatrix("dehomogenize_data.csv");
h = data(1);
nCells = data(2);
radius = data(3:nCells+2) .* 1e-3;
center = reshape(data(nCells+3:end),[],2);


%% Draw Cells in 2D and 3D
figure(1); clf(1);
for i = 1:length(center)
    drawCell(radius(i),center(i,:),h);
end
axis equal;
xlim([0,0.1]); ylim([0,0.05]);
title("Basic Dehomogenization in 2D");

figure(2);
nLayers = 5;
for j = 0:nLayers-1
    for i = 1:length(center)
        drawCell3d(radius(i), [center(i,1),j*h,center(i,2)], h);
    end
end
axis equal;

function drawCell(r, c, h)
    A = [c(1), c(2)+0.5*h];
    B = [c(1)-0.5*h, c(2)];
    C = [c(1), c(2)-0.5*h];
    D = [c(1)+0.5*h, c(2)];
    dx = [sqrt(2*r^2),0];
    dy = [0,sqrt(2*r^2)];

    points = [A-dx;
        A+dx;
        D+dy;
        D-dy;
        C+dx;
        C-dx;
        B-dy;
        B+dy;
        A-dx];
    pointsIn = [A-dy;
        D-dx;
        C+dy;
        B+dx;
        A-dy];

    plot(points(:,1),points(:,2),"Color",[0 0.4470 0.7410],"LineWidth",2.0); hold on;
    plot(pointsIn(:,1),pointsIn(:,2),"Color",[0 0.4470 0.7410],"LineWidth",2.0);
end

function drawCell3d(r,c,h)
    cellCross = sqrt(2 * (h^2));
    cellCrossThrough = sqrt(h^2 + cellCross^2);
    strutLength = cellCrossThrough / 2;
    strutAngle = atand(cellCross/h);
    alpha = (0:pi/10:2*pi)';
    circlePoints = [0,0;r * sin(alpha), r * cos(alpha)];
    cylinderPoints = [[circlePoints, zeros(22,1)]; [circlePoints, strutLength * ones(22,1)]];

    conn = [2,3,24;3,24,25];
    ch = (0:19)' .* ones(20,3);
    connfull = [conn(1,:) + ch;conn(2,:) + ch];
    connfull = repmat(connfull,8,1);
    connfull(41:80,:) = connfull(1:40,:) + 44;
    connfull(81:120,:) = connfull(41:80,:) + 44;
    connfull(121:160,:) = connfull(81:120,:) + 44;
    connfull(161:200,:) = connfull(121:160,:) + 44;
    connfull(201:240,:) = connfull(161:200,:) + 44;
    connfull(241:280,:) = connfull(201:240,:) + 44;
    connfull(281:320,:) = connfull(241:280,:) + 44;

    strut1 = rotateZ(rotateY(cylinderPoints, strutAngle),45);
    strut2 = rotateZ(rotateY(cylinderPoints, strutAngle),45+90);
    strut3 = rotateZ(rotateY(cylinderPoints, strutAngle),45+180);
    strut4 = rotateZ(rotateY(cylinderPoints, strutAngle),-45);
    struts = [strut1;strut2;strut3;strut4];
    struts5 = rotateY(struts, 180);
    struts = [struts;struts5];
    struts = struts + ones(length(struts),3) .* c;

    trisurf(connfull, struts(:,1),struts(:,2),struts(:,3),'FaceColor','blue','LineStyle','none'); hold on;
end

function p = rotateY(points, angle)
    nPoints = size(points,1);
    rotY = [cosd(angle), 0, -sind(angle);
        0, 1, 0;
        sind(angle), 0, cosd(angle);];

    p = zeros(size(points));
    for i = 1:nPoints
        p(i,:) = (rotY * points(i,:)')';
    end
end

function p = rotateZ(points, angle)
    nPoints = size(points,1);
    rotY = [cosd(angle), -sind(angle), 0;
        sind(angle), cosd(angle), 0;
        0, 0, 1];

    p = zeros(size(points));
    for i = 1:nPoints
        p(i,:) = (rotY * points(i,:)')';
    end
end



