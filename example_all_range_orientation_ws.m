% Author:     Changmeng Hou(Harbin Engineering University)

% All Range Orientation Workspace.

clc
clear
tic

RDL = [5 3 .5 .3 7.2 10];
[B,P,H] = rdl2bph(RDL);
limbPara = RDL(5:6);
jointCons = [0,pi/4];

euler = [0 0 0];
almax = pi/18;
bemax = pi/18;
gamax = pi/18;
% disStep = 0.02;
disStep = 0.04;
eu_i = 1;
euler(eu_i,:) = [0,0,0];
for al = -almax:disStep:almax
    for be = -bemax:disStep:bemax
        for ga = -gamax:disStep:gamax
            eu_i = eu_i + 1;
            euler(eu_i,:) = [al,be,ga];
        end
    end
end


zmin = H(1);
zmax = H(2);
% numz = 50;
% numtheta = 90;
numz = 40;
numtheta = 60;
z = linspace(zmin,zmax,numz);
theta = linspace(0,2*pi,numtheta);

deltaz = (zmax-zmin)/numz;
deltatheta = 2*pi/numtheta;
deltarho = 0.02;

% Z = zeros(numz,numtheta);
Rhotemp = zeros(numz,numtheta);
% Theta = zeros(numz,numtheta);
% XC = zeros(numz,numtheta);
% YC = zeros(numz,numtheta); 
% ZC = zeros(numz,numtheta);
tic

subvi = 0;
subv = zeros(numz*numtheta,1);
for i =1:numz
    for j = 1:numtheta
        rho = 0;
        while (rho >= 0)
            [xc,yc,zc] = pol2cart(theta(j),rho,z(i));
            transTemp = [xc,yc,zc];
            for ei = 1:eu_i
                eulerTemp = euler(ei,:);
                posTemp = [transTemp,eulerTemp];
                [inWorkspace,isSingular,inBound] = pose_state(posTemp,jointCons,limbPara,B,P);
                if ~inWorkspace
                    break
                end
                
            end
            if ~inWorkspace
                Rhotemp(i,j) = rho;
                subvi = subvi + 1;
                subv(subvi) = 1/2*(rho^2)*deltaz*deltatheta;
                break
            end
            rho = rho + deltarho;

        end
    end
toc
end


[row] = find(all(Rhotemp ~= 0,2));
if (max(row) ~= numz)&&(min(row) ~= 1)
    row = [min(row)-1;row;max(row)+1];
elseif (max(row) == numz)&&(min(row) ~= 1)
    row = [min(row)-1;row];
elseif (max(row) ~= numz)&&(min(row) == 1)
    row = [row;max(row)+1];
else
    
end

Rho = Rhotemp(row,:);
z = z(row);
[Theta,Z] = meshgrid(theta,z);
[XC,YC,ZC] = pol2cart(Theta,Rho,Z);
mesh(XC,YC,ZC);
workspaceV = sum(subv)
toc



