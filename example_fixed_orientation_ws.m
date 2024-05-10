% Author:     Changmeng Hou(Harbin Engineering University)

% Fixed Orientation - Location Workspace.

clc
clear
tic

RDL = [5 3 .5 .3 7.2 10];
[B,P,H] = rdl2bph(RDL);
limbPara = RDL(5:6);
jointCons = [0,pi/4];


zmin = H(1);
zmax = H(2);
zRange = [zmin,zmax];
euler = [0 0 0];
numz = round((zRange(2)-zRange(1))/0.1);
numTheta = 90;
z = linspace(zRange(1),zRange(2),numz);

theta = linspace(0,2*pi,numTheta);

% theta = distributed(theta);
deltaz = (zRange(2)-zRange(1))/numz;
deltaTheta = 2*pi/numTheta;
deltaRho = 0.02;

rhoTemp = zeros(numz,numTheta);
subVoli = 0;
subVol = zeros(numz*numTheta,1);
V = 0;

for i =1:numz
    for j = 1:numTheta
        rho = 0;
        inWorkspace = 1;
        while (rho >= 0&& inWorkspace)
            
            [xc,yc,zc] = pol2cart(theta(j),rho,z(i));
            transTemp = [xc,yc,zc];
            eulerTemp = [0 0 0];
            poseTemp = [transTemp,eulerTemp];
            [inWorkspace,isSingular,inBound] = pose_state(poseTemp,jointCons,limbPara,B,P);
            rho = rho + deltaRho;
        end
        rhoTemp(i,j) = rho ;
        subV = 1/2*(rho^2)*deltaz*deltaTheta;
        V =V+subV;
    end
end


[row] = find(all(rhoTemp ~= 0,2));
if (max(row) ~= numz)&&(min(row) ~= 1)
    row = [min(row)-1;row;max(row)+1];
elseif (max(row) == numz)&&(min(row) ~= 1)
    row = [min(row)-1;row];
elseif (max(row) ~= numz)&&(min(row) == 1)
    row = [row;max(row)+1];
else
    
end

RHO = rhoTemp(row,:);
z = z(row);
[THETA,Z] = meshgrid(theta,z);
[XC,YC,ZC] = pol2cart(THETA,RHO,Z);
mesh(XC,YC,ZC)


toc

