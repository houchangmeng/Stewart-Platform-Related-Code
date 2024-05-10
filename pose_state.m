% Author:     Changmeng Hou(Harbin Engineering University)

% Return current platform state(1. is in workspace ? 2. is singular? 3. is in boundary?).

function [inWorkspace,isSingular,inBound] = pose_state(pose,jointCons,limbCons,B,P)


[nowThetaPi,nowThetaBi,limbLength,nowDMin,Jacobi] = constrained_state(pose,B,P);

% is singularity ?
isSingular = 0;
if det(Jacobi) == 0
    isSingular = 1;
end

thetaPiInterval = [min(nowThetaPi),max(nowThetaPi)];
thetaBiInterval = [min(nowThetaBi),max(nowThetaBi)];
limbInterval = [min(limbLength),max(limbLength)];
dMin = min(nowDMin);

thetaPiCons = jointCons;
thetaBiCons = jointCons;

overPi = 0;
overBi = 0;
overLimb = 0;
overD = 0;
% is satisfy the joints constraints?
if (thetaPiInterval(1)<thetaPiCons(1))||(thetaPiInterval(2)>thetaPiCons(2))
    overPi = 1;
end
if (thetaBiInterval(1)<thetaBiCons(1))||(thetaBiInterval(2)>thetaBiCons(2))
    overBi = 1;
end
% is satisfy the limb length constraints?
if (limbInterval(1)<limbCons(1))||(limbInterval(2)>limbCons(2))
    overLimb = 1;
end
% is satisfy the limb length constraints?
if dMin < 0.1
    overD = 1;
end

inWorkspace = 1;
if overPi||overBi||overLimb||overD
    inWorkspace = 0;
end

boundPi = 0;
boundBi = 0;
boundLimb = 0;
boundD = 0;
% in workspace bound?
if inWorkspace
    if ((thetaPiInterval(1)-thetaPiCons(1)) < 0.05)||((thetaBiCons(2)-thetaPiInterval(2)) < 0.05)
        boundPi = 1;
    end
    if ((thetaBiInterval(1)-thetaBiCons(1)) < 0.05)||((thetaBiCons(2)-thetaBiInterval(2)) < 0.05)
        boundBi = 1;
    end
    if (limbInterval(1) - limbCons(1) < 0.2)||((limbCons(2)-limbInterval(2)) < 0.2)
        boundLimb = 1;
    end
        
        boundD = 0;
end
% in bound?
inBound = 0;
if boundPi||boundBi||boundLimb||boundD
    inBound = 1;
end


end

% Return current joint location , limb length and jacobian.
function [nowThetaPi,nowThetaBi,nowLimbLength,nowDMin,Jacobi] = constrained_state(pose,B,P)

%   theta_pi : angles of move platform joints.（deg）
%   theta_bi : angles of fixed platform joints.（deg）  
%   npi : norm vector of ith joint in move platform.
%   nbi : norm vector of ith joint in fixed platform
%   N : public norm vector of adjacent limb vector
%   M(1:3,i): mi
%   M(4:6,i): mi+1
%   CB : vector between public norm vector and fixed platform of adjacent limb
%   cb : distance between public norm vector and fixed platform of adjacent limb
%   d : the shortest distance of of  adjacent limbs.
%   dmin : the shortest distance of all limbs;
%   mb(1,i) : projection of limb(i+1) in limb(i)
%   mb(2,i) : projection of limb(i) in limb(i+1)
 
% position inverse solution
[limbLength,Q,L,R_PB,Jacobi] = IK(pose,B,P);
nowLimbLength = limbLength;
% Extend to [ i ; i+1 ]
limb_1 = nowLimbLength;        % i
limb_2 = [nowLimbLength(2:6),nowLimbLength(1)];
                            % i+1
limb_ex = [limb_1;
          limb_2];
B_ex = [B;
        B(:,2:6),B(:,1)];
Q_ex = [Q;                  % i
        Q(:,2:6),Q(:,1)];   % i+1
L_ex = [L;
        L(:,2:6),L(:,1)];

% Initial state of ball joints
npi = [[0 0 0 0 0 0];[0 0 0 0 0 0];[-1 -1 -1 -1 -1 -1]];
nbi = [[0 0 0 0 0 0];[0 0 0 0 0 0];[ 1  1  1  1  1  1]];

% % Angle of ball joints
if ~isa(pose,'sym')
    nowThetaPi = zeros(1,6);
    nowThetaBi = zeros(1,6);
end

for i = 1:6
    nowThetaPi(i) = vector_angle(-L_ex(1:3,i),R_PB*npi(:,i)); 
    nowThetaBi(i) = vector_angle(L_ex(1:3,i),nbi(:,i));
end


% % mechanical interference 
nowDMin = ones(1,6);

% Bextend = [B;B(:,1)];
% Qextend = [Q;Q(:,1)];
% for iDist = 1:6
%     p1 = Bextend(:,iDist);
%     p2 = Qextend(:,iDist);
%     p3 = Bextend(:,iDist+1);
%     p4 = Qextend(:,iDist+1);
%     currentLimbDistance(iDist) = DistBetween2Segment(p1, p2, p3, p4);
% end
% nowDMin = currentLimbDistance;


% Method in <<高等空间机构学>>
% N  = zeros(6,6);         
% m  = zeros(6,6);         
% CB = zeros(6,6);       
% mb = zeros(2,6);                  
% d  = zeros(1,6);
% cb = zeros(2,6);
% b  = zeros(6,6);
% bp = zeros(1,6);


% % 1~3 ni,mi,ci, 4~6 ni+1 mi+1 ci+1
% 
% for i = 1:6
%     N(1:3,i) = cross(L_ex(1:3,i),L_ex(4:6,i))/norm(cross(L_ex(1:3,i),L_ex(4:6,i))); 
%     N(4:6,i) = cross(L_ex(4:6,i),L_ex(1:3,i))/norm(cross(L_ex(4:6,i),L_ex(1:3,i)));
%     d(i) = abs(dot(N(1:3,i),B_ex(4:6,i)-B_ex(1:3,i)));
%     m(1:3,i) = cross(N(1:3,i),L_ex(4:6,i)); % mi
%     m(4:6,i) = cross(N(4:6,i),L_ex(1:3,i)); % mi+1
%     
%     bp(i) = (norm(Q_ex(1:3,i)-Q_ex(4:6,i))) < (norm(B_ex(1:3,i)-B_ex(4:6,i)));
%     if bp(i)
% 
%         b(1:3,i) = B_ex(1:3,i);
%         b(4:6,i) = B_ex(4:6,i);
%     else
% 
%         b(1:3,i) = Q_ex(1:3,i);
%         b(4:6,i) = Q_ex(4:6,i);
%     end
%     CB(1:3,i) = abs(dot(b(4:6,i)-b(1:3,i),m(1:3,i))/dot(L_ex(1:3,i),m(1:3,i)))*L_ex(1:3,i);%ci-bi
%     CB(4:6,i) = abs(dot(b(1:3,i)-b(4:6,i),m(4:6,i))/dot(L_ex(4:6,i),m(4:6,i)))*L_ex(4:6,i);%ci+1-bi+1
%     
%     cb(1,i) = norm(CB(1:3,i));  % cb_i
%     cb(2,i) = norm(CB(4:6,i));  % cb_i+1
%     mb(1,i) = abs(dot(L_ex(1:3,i),Q_ex(4:6,i)-B_ex(1:3,i)));    % mb_i
%     mb(2,i) = abs(dot(L_ex(4:6,i),Q_ex(1:3,i)-B_ex(4:6,i)));    % mb_i+1
%     
%     % % The shortest distance of adjacent two limb.
%     
%     if (cb(1,i)<limb_ex(1,i))&&(cb(2,i)<limb_ex(2,i)) 
%         % case 1
%         
%         dmin(i) = d(i);
%     elseif (cb(1,i)>limb_ex(1,i))&&(cb(2,i)<limb_ex(2,i))
%         % case 2
%         
%         % Page-163  (6-92)
%         dmin(i) = norm(cross(Q_ex(1:3,i)-B_ex(4:6,i),L_ex(4:6,i)))/norm(L_ex(4:6,i));
%     elseif (cb(1,i)<limb_ex(1,i))&&(cb(2,i)>limb_ex(2,i)) 
%         
%         % Page-163  (6-93)
%         dmin(i) = norm(cross(Q_ex(4:6,i)-B_ex(1:3,i),L_ex(1:3,i)))/norm(L_ex(1:3,i));
%     else 
%         % case 3
%         if (mb(1,i)>limb_ex(1,i))&&(mb(2,i)<limb_ex(2,i))
%             % Page-163  ①
%             % Page-163  (6-92)
%             dmin(i) = norm(cross(Q_ex(1:3,i)-B_ex(4:6,i),L_ex(4:6,i)))/norm(L_ex(4:6,i));
%         elseif (mb(1,i)<limb_ex(1,i))&&(mb(2,i)>limb_ex(2,i))
%             % Page-163  ②
%             % Page-163  (6-93)
%             dmin(i) = norm(cross(Q_ex(4:6,i)-B_ex(1:3,i),L_ex(1:3,i)))/norm(L_ex(1:3,i));
%         elseif  (mb(1,i)>limb_ex(1,i))&&(mb(2,i)>limb_ex(2,i))
%             % Page-163  ③
%             if bp(i)
%                 dmin(i) = norm(Q_ex(4:6,i)-Q_ex(1:3,i));
%             else
%                 dmin(i) = norm(B_ex(4:6,i)-B_ex(1:3,i));
%             end
%         else
% 
%         end
%     end
% end


end
function angleRad = vector_angle(v1,v2)
% vector_angle calculate the angle of 2 vectors(rad)   

angleRad = acos(dot(v1,v2) / ( norm(v1)*norm(v2)) );

end

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



