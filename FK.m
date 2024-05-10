% Author:     Changmeng Hou(Harbin Engineering University)

% Forward Kinematics.

function [pose,isSingular] = FK(lmbLength,B,P)

if (~isequal(size(P),[3,6]))||(~isequal(size(B),[3,6]))
    error('Check the paramaters P or B!')
end
if (~isequal(length(lmbLength),6))
    error('Check the paramaters limblength')
end

% initial position
initPose = [0,0,mean(lmbLength),0,0,0].';

finalLmbL = lmbLength;

poseTemp = initPose(:);
[nowLmbL,~,~,~,J] = IK(poseTemp,B,P);
deltaLmb = finalLmbL(:)-nowLmbL(:);
Jmodi = getJmodi(poseTemp,J);
pose = poseTemp+Jmodi\deltaLmb;

iteraIdx = 0;

while (norm(pose - poseTemp)>1e-3) && (iteraIdx < 50)
    
    poseTemp = pose;
    [nowLmbL,~,~,~,J] = IK(poseTemp,B,P);
    % delta limb
    deltaLmb = finalLmbL(:)-nowLmbL(:);
    Jmodi = getJmodi(poseTemp,J);
    % singular judge
    poseVel = Jmodi\deltaLmb;  
    isSingular = norm(poseVel,'inf')>(100*mean(finalLmbL)); 
    if isSingular
        break
    end
    pose = poseTemp+poseVel;
    iteraIdx = iteraIdx +1;
end

function Jmodi = getJmodi(pose,J)

be = pose(5);
ga = pose(6);
sb = sin(be);
sg = sin(ga);
cb = cos(be);
cg = cos(ga);

U = [cg*cb  -sg  0;
     sg*cb   cg  0;
     -sb     0   1];
Jmodi = J*[ eye(3),  zeros(3);
            zeros(3),   U    ]; 
end
end