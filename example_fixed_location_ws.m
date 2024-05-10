% Author:     Changmeng Hou(Harbin Engineering University)

% Fixed Location - Orientation Workspace.

clc
clear
tic

RDL = [5 3 .5 .3 7.2 10];
[B,P,H] = rdl2bph(RDL);
limbPara = RDL(5:6);
jointCons = [0,pi/4];


gammarange = [-pi/2,pi/2];

pose = FK(repmat(mean(limbPara),[1,6]),B,P);
transTemp = pose(1:3).';

numGamma = round((pi/0.05));
numTheta = 90;
gamma = linspace(-pi/2,pi/2,numGamma);
theta = linspace(0,2*pi,numTheta);

deltaGamma = pi/numGamma;
deltaTheta = 2*pi/numTheta;
deltaRho = 0.005;

rhoTemp = zeros(numGamma,numTheta);
subVoli = 0;
subVol = zeros(numGamma*numTheta,1);
for i =1:numGamma
    for j = 1:numTheta
        rho = 0;
        while (rho >= 0)
            
            [alc,bec,gac] = pol2cart(theta(j),rho,gamma(i));
%             transTemp = [0 0 .5*(H(1)+H(2))];
            eulerTemp = [alc,bec,gac];
            poseTemp = [transTemp,eulerTemp];
            [inWorkspace,isSingular,inBound] = pose_state(poseTemp,jointCons,limbPara,B,P);
            if inWorkspace
                rho = rho + deltaRho;
            else
                rhoTemp(i,j) = rho ;
                subVoli = subVoli + 1;
                subVol(subVoli) = 1/2*(rho^2)*deltaGamma*deltaTheta;
                break
            end
        end
    end
end


[row] = find(all(rhoTemp ~= 0,2));
if (max(row) ~= numGamma)&&(min(row) ~= 1)
    row = [min(row)-1;row;max(row)+1];
elseif (max(row) == numGamma)&&(min(row) ~= 1)
    row = [min(row)-1;row];
elseif (max(row) ~= numGamma)&&(min(row) == 1)
    row = [row;max(row)+1];
else
    
end

RHO = rhoTemp(row,:);
gamma = gamma(row);
[THETA,Z] = meshgrid(theta,gamma);
[XC,YC,ZC] = pol2cart(THETA,RHO,Z);

mesh(XC,YC,ZC)
% COWVolume = sum(subVol)
% XCdeg = XC.*180./pi;
% YCdeg = YC.*180./pi;
% ZCdeg = ZC.*180./pi;
% mesh(XCdeg,YCdeg,ZCdeg)

toc


