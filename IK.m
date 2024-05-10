% Author:     Changmeng Hou(Harbin Engineering University)

% Inverse Kinematics.

function [nowLmbL,Q,L,R_PB,Jac] = IK(pose,B,P)

if ~( (isequal(length(pose),6))&&(isequal(size(B),[3,6]))&&(isequal(size(P),[3,6])) )
    error('Check your input!')
end
    
Trans = pose(1:3);
Euler = pose(4:6);

% % calculate the output

T = Trans(:);
R_PB = rotz(Euler(3))*roty(Euler(2))*rotx(Euler(1)); % rotation (zyx)
P_B = R_PB*P;
Q = T+P_B;
L = Q-B;

if ~isa(pose,'sym')
    nowLmbL = zeros(1,6);
    s = zeros(3,6);
    Jac = zeros(6,6);
end

for i = 1:6    
    nowLmbL(i) = norm(L(:,i));
    s(:,i) = L(:,i)/nowLmbL(i);
    Jac(i,:) = [s(:,i).' cross(P_B(:,i),s(:,i)).'];
end

detJ = det(Jac); 

end

