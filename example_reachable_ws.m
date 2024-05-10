% Author:     Changmeng Hou(Harbin Engineering University)

% Reachable Workspace in Orientation Range.

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
rhoMax = 4;

alMax = pi/12;
beMax = pi/12;
gaMax = pi/12;
step = 8;
discreteStep = alMax/step;
eu_i = 0;
% initEulerGroup(eu_i,:) = [0 0 0];
for al = -alMax:discreteStep:alMax
    for be = -beMax:discreteStep:beMax
        for ga = -gaMax:discreteStep:gaMax
            eu_i = eu_i + 1;
            initEulerGroup(eu_i,:) = [al,be,ga];
        end
    end
end
NewEulerGroup = [[0 0 0];initEulerGroup];
numz = 30;
numTheta = 30;
numDelta = 30;
z = linspace(zRange(1),zRange(2),numz);
theta = linspace(0,2*pi,numTheta);

deltaz = (zRange(2)-zRange(1))/numz;
deltaTheta = 2*pi/numTheta;
deltaRho = rhoMax/numDelta;

rhoTemp = zeros(numz,numTheta);
subVoli = 0;
subVol = zeros(numz*numTheta,1);
tic
zSearchPoint = 0;
V = 0;
for i =1:numz
    for j = 1:numTheta
        rho = 0;
        while (rho < rhoMax)

            [xc,yc,zc] = pol2cart(theta(j),rho,z(i));
            transTemp = [xc,yc,zc];
%             NewEulerGroup = getEulerGroup(theta(j),zSearchPoint,initEulerGroup);

                for ei = 1:eu_i
                    eulerTemp = NewEulerGroup(ei,:);
                    posTemp = [transTemp,eulerTemp];
                    [inWorkspace,isSingular,inBound] = pose_state(posTemp,jointCons,limbPara,B,P);
                    if inWorkspace
                        break
                    end
                end
                    if ~inWorkspace
                        rhoTemp(i,j) = rho ;
                        subV = 1/2*(rho^2)*deltaz*deltaTheta;
                        V =V+subV;
%                         subVoli = subVoli + 1;
%                         subVol(subVoli) = 1/2*(rho^2)*deltaz*deltaTheta;
                        break
                    end
            rho = rho + deltaRho;
        end
    end
%     i
    toc
end
toc

[row] = find(all(rhoTemp ~= 0,2));
if (max(row) ~= numz)&&(min(row) ~= 1)
    row = [min(row)-1;row;max(row)+1];
elseif (max(row) == numz)&&(min(row) ~= 1)
    row = [min(row)-1;row];
elseif (max(row) ~= numz)&&(min(row) == 1)
    row = [row;max(row)+1];
else
    % do nothing
end

RHO = rhoTemp(row,:);
z = z(row);
[THETA,Z] = meshgrid(theta,z);
[XC,YC,ZC] = pol2cart(THETA,RHO,Z);
mesh(XC,YC,ZC)
V



