% Author:     Changmeng Hou(Harbin Engineering University)
% Reference : Lung-Wen Tsai
%             Solving the Inverse Dynamics of a Stewart-Gough Manipulator 
%             by the Principle of Virtual Work

clc
clear
tic
                        % % Initial condition % %

P = [0.17 0.595 -0.4; -0.6 0.15 -0.4; -0.6 -0.15 -0.4; 0.17 -0.595 -0.4;0.43 -0.445 -0.4;0.43 0.445 -0.4]';
B = [-2.12 1.374 0; -2.38 1.224 0; -2.38 -1.224 0;-2.12 -1.374 0;0 -0.15 0; 0 0.15 0]';
mp = 1.5;
m1 = 0.1;
m2 = 0.1;
g = [0 0 -9.80665]';
I_P = diag([0.08 0.08 0.08]);
I1_i = diag([1/160 1/160 1/160]);
I2_i = diag([1/160 1/160 1/160]);
e1 = 0.5;
e2 = 0.5;
initL = [0 0 0]; 

% input forces exerted at the center of mass of the platform
fe = [0 0 0]';
ne = [0 0 0]';

syms t;
assume(t,'Real')
assumeAlso(t>0)

w = 3; % unit:rad/s
                        
% % Motion pattern 1  
 Trans = [-1.5+0.2*sin(w*t),0.2*sin(w*t),1+0.2*sin(w*t)].';
 Euler = [sym(0),sym(0),sym(0)].';
 
% % Motion pattern 2
%  Trans = [-1.5,0,1]';
%  Euler = [sym(0),sym(0),pi/6*sin(w*t)]';
 
% 
% % Motion pattern 3 
%  Trans = [-1.5+0.2*sin(w*t),0,1]';
%  Euler = [sym(0),sym(0),sym(0)]';

% % Motion pattern 4 
%  Trans = [-1.5,0,1].';
%  Euler = [0,0,0].';

T = Trans;              % position of the platform's mass center
vt = diff(T,t);         % velocity of the platform's mass center
at = diff(vt,t);        % acceleration of the platform's mass center

R_PB = rotz(Euler(3))*roty(Euler(2))*rotx(Euler(1));
                        % Rotation matrix of frame platform relative to frame base(zyx)
P_B = R_PB*P;           % position of (i)st upper ball joint in the fixed frame
L = T + P_B - B;        % postion inverse solution

                        % % % Allocate space % % %
                        
d = cell(1,6);          % limb length
s = cell(1,6);          % unit vector of limb from fixed frame to moving frame
R_iB = cell(1,6);       % rotation matrix of limb frame i relative to fixed frame
R_Bi = cell(1,6);       % rotation matrix of fixed frame relative to limb frame i

                        % % velocity % %
                    
vp_B = cell(1,6);       % velocity of upper joint Pi
vp_i = cell(1,6);       % velocity of upper joint Pi in the ith limb frame
vL_i_z = cell(1,6);     % limb's velocity in the limb frame(z axis)
omegaL_i = cell(1,6);   % angular velocity of the ith limb in the ith limb frame
v1_i = cell(1,6);       % velocity of the under limb's in the ith limb frame 
v2_i = cell(1,6);       % velocity of the upper limb's in the ith limb frame

                        % % acceleration % %
                    
[omegat,alphat] = euler2omealpvec(Euler); 
                        % Euler angular velocity is transformed into angular velocity vector and angular acceleration vector
                        % omega_p:  angular velocity of the mass center of moving platform
                        % alpha_p:  angular acceleration of the mass center of moving platform
                        
ap_B = cell(1,6);       % acceleration of the upper joint
ap_i = cell(1,6);       % acceleration of the upper joint in the limb frame
aL_i_z = cell(1,6);     % acceleration of the ith limb in the limb frame(z-axis)
alphaL_i = cell(1,6);   % angular acceleration of the ith limb in the limb frame
a1_i = cell(1,6);       % acceleration of the under limb's in the limb frame
a2_i = cell(1,6);       % acceleration of the upper limb's in the limb frame

                        % % Jacobi % %                
                    
Jb= cell(1,6);          % Jacobi of the ith joint 
Jb_i= cell(1,6);        % Jacobi of the ith limb

J1_i = cell(1,6);       % Jacobi of the ith under limb 
J2_i = cell(1,6);       % Jacobi of the ith upper limb 

JP = sym('JL',[6 6]);                      
Jx = sym('Jx',[6 6]);   % Jacobi of the ith limbs(limbframe x axis)
Jy = sym('Jy',[6 6]);   % Jacobi of the ith limbs(limbframe y axis)
 
                        % % Force and Movement % %
               
f1 = cell(1,6);         % resulant of applied and inertia forces exerted at the ith under limb
n1 = cell(1,6);         % resulant of applied and inertia moments exerted at the  ith under limb
F1_i = cell(1,6);       % wrenches F1 = [f1;n1]
                        
f2 = cell(1,6);         % resulant of applied and inertia forces exerted at the ith upper limb
n2 = cell(1,6);         % resulant of applied and inertia moments exerted at the ith upper limb
F2_i = cell(1,6);       % wrenches F2 = [f2;n2]

I_B = R_PB*I_P*R_PB';
fp = fe+mp*g-mp*at;     % resulant of applied and inertia forces exerted at the mass center of platform
np = ne-I_B*alphat-cross(omegat,(I_B*omegat));
                        % resulant of applied and inertia moments exerted at the mass center of platform
Fp = [  fp  ;
        np  ];    

Fx = sym('Fx',[6 1]);   % the forces in the x-axis of all limbs in the limb frame
Fy = sym('Fy',[6 1]);   % the forces in the y-axis of all limbs in the limb frame
Fz = sym('Fz',[6 1]);   % the forces in the z-axis of all limbs in the limb frame
                        
for i = 1:6
    
    % limb length
    d{i} = norm(L(:,i));
    % unit vector of limb
    s{i} = L(:,i)/d{i};

    % rotation matrix of limb frame i relative to fixed frame
    ct = s{i}(3);
    st = sqrt(s{i}(1)^2+s{i}(2)^2);
    sf = s{i}(2)/st;
    cf = s{i}(1)/st; 
    R_iB{i} = [
               -cf*ct -sf  cf*st;
                sf*ct  cf  sf*st;
               -st     0    ct  ;
              ];
    R_Bi{i} = R_iB{i}.';
    
                        % % velocity % %
                        
     % velocity of the ith upper joint
     vp_B{i} = vt+cross(omegat,P_B(:,i));
     vp_i{i} = R_Bi{i}*vp_B{i};
     % velocity of the ith limb in the ith limb frame (z-axis)
     vL_i_z{i} = vp_i{i}(3);
     % angular velocity of the ith limb in the ith limb frame
     omegaL_i{i} = 1/d{i}*(cross([0 0 1]',vp_i{i}));
     % velocity of the under limb's in the ith limb frame 
     v1_i{i} = e1*cross(omegaL_i{i},[0,0,1]');
     % velocity of the upper limb's in the ith limb frame 
     v2_i{i} = (d{i} - e2)*cross(omegaL_i{i},[0,0,1]')+vL_i_z{i}*[0 0 1]';

                        % % acceleration % %
                        
     % acceleration of the ith upper joint
     ap_B{i} = at + cross(alphat,P_B(:,i))+cross(omegat,cross(omegat,P_B(:,i)));
     ap_i{i} = R_Bi{i}*ap_B{i};
     % acceleration of the ith limb in the ith limb frame (z-axis)
     aL_i_z{i} = diff(vL_i_z{i},t) + d{i}*dot(omegaL_i{i},omegaL_i{i});
     % angular acceleration of the ith limb in the limb frame
     alphaL_i{i} = 1/d{i}*cross([0 0 1]',ap_i{i})...
                         -2*vL_i_z{i}/d{i}*omegaL_i{i};               
     % acceleration of the ith under limb's in the limb frame           
     a1_i{i} = e1*cross(alphaL_i{i},[0 0 1]') + e1*cross(omegaL_i{i},cross(omegaL_i{i},[0 0 1]'));
     % acceleration of the ith upper limb's in the limb frame 
     a2_i{i} = aL_i_z{i}*[0 0 1]'+(d{i} - e2)*cross(alphaL_i{i},[0 0 1]')+...
                    (d{i}-e2)*cross(omegaL_i{i},cross(omegaL_i{i},[0 0 1]'))+...
                     2*vL_i_z{i}*cross(omegaL_i{i},[0 0 1]');
                 
                        % % force moments wrenches %%
                        
     % resulant of applied and inertia forces exerted at the ith under limb
     f1{i} = m1*R_Bi{i}*g - m1*a1_i{i};
     % % moments 
     n1{i} = -I1_i*alphaL_i{i} - cross(omegaL_i{i},(I1_i*omegaL_i{i}));
     % % wrenches 
     F1_i{i} = [
                 f1{i};
                 n1{i}
                ];
     % resulant of applied and inertia forces exerted at the ith upper limb
     f2{i} = m2*R_Bi{i}*g - m2*a2_i{i};
     % % moments
     n2{i} = -I2_i*alphaL_i{i} - cross(omegaL_i{i},(I2_i*omegaL_i{i}));
     % % wrenches
     F2_i{i} = [
                 f2{i};
                 n2{i}
                ]; 
     % the forces in the x-axis of all limbs in the limb frame
     Fx(i,1) = 1/d{i}*(e1*f1{i}(1)+(d{i}-e2)*f2{i}(1)+n1{i}(2)+n2{i}(2));
     % % y-axis
     Fy(i,1) = 1/d{i}*(e1*f1{i}(2)+(d{i}-e2)*f2{i}(2)+n1{i}(1)+n2{i}(1));
     % % z-axis
     Fz(i,1) = f2{i}(3);
     
                        % % Jacobi % %           

     % Jacobi of the ith upper joint
     Jb{i} = [eye(3),-skew(P_B(:,i))];
     % Jacobi of the ith limb(in the limb frame)
     Jb_i{i} = R_Bi{i}*Jb{i}; 
     % Jacobi of manipulator
     JP(i,:) = Jb_i{i}(3,:);
     % Jacobi of the ith limbs(limb frame x-axis)
     Jx(i,:) = Jb_i{i}(1,:);
     % % y-axis
     Jy(i,:) = Jb_i{i}(2,:);
    
end
toc

                        % % plot figure % %
                        
f = figure('Position',[560 75  475  680]);
%
timeInterval = [0,4*pi/w];
% 
subplot(4,1,1);
fplot(d,timeInterval);
legend('d_1','d_2','d_3','d_4','d_5','d_6');
title('Leg length')
xlabel('t(s)')
ylabel('length(m)')

subplot(4,1,2)
fplot(vL_i_z,timeInterval);
legend('v_1','v_2','v_3','v_4','v_5','v_6');
title('Velocity')
xlabel('t(s)')
ylabel('velocity(m/s)')

subplot(4,1,3)
fplot(aL_i_z,timeInterval);
legend('a_1','a_2','a_3','a_4','a_5','a_6');
title('Accleration')
xlabel('t(s)')
ylabel('accleration(m/s^2)')

toc
%   the first dynamics inverse solution

% Symbol calculation takes too long, so convert it to numerical value to speed up the calculation
% tau1 = -Fz - (Jz.')\(Fp+(Jx.')*Fx+(Jy.')*Fy); % Fmula

n = 20; 
T = linspace(0,timeInterval(2),n); 
tau1 = [];
for tt = T
    
    FxTemp = double(subs(Fx,t,tt));
    FyTemp = double(subs(Fy,t,tt));
    FzTemp = double(subs(Fz,t,tt));
    FpTemp = double(subs(Fp,t,tt));
    
    JxTemp = double(subs(Jx,t,tt));
    JyTemp = double(subs(Jy,t,tt));
    JpTemp = double(subs(JP,t,tt));
    
    % Calculate the force at the current moment using the same formula as noted above  
    tauTemp = -FzTemp - (JpTemp.')\(FpTemp+(JxTemp.')*FxTemp+(JyTemp.')*FyTemp);
    tau1 = [tau1,tauTemp];
end
toc

subplot(4,1,4)
plot(T,tau1.');
legend('f_1','f_2','f_3','f_4','f_5','f_6')
title('Dynamics')
xlabel('t(s)')
ylabel('Force(N)')

function [omega_p,alpha_p] = euler2omealpvec(euler)
% This function convert Euler  to omega and alpha vectors
% Syntax :
%           [omega_p,alpha_p] = euler2omealpvec(Euler)
% Input :
%           Euler - the euler angle
% Output :
%           omega_p - the angle velocity vector in this Euler
%           alpha_p - the angle accleration vector in this Euler   
% Reference ： Liu's Book

syms t;
assume(t,'Real')
assumeAlso(t>0)

    alpha = euler(1);
    beta  = euler(2);
    gamma = euler(3);
    
    cb = cos(beta);
    cg = cos(gamma);

    sb = sin(beta);
    sg = sin(gamma);

    da = diff(alpha,t);
    db = diff(beta,t);
    dg = diff(gamma,t);

    dda = diff(da,t);
    ddb = diff(db,t);
    ddg = diff(dg,t);
    

    U = [
         cg*cb  -sg  0;
         sg*cb   cg  0;
         -sb      0  1
         ];
    omega_p = U*[da db dg]';
    W = [
         -dg*sg*cb-db*cg*sb  -dg*cg  0;
          dg*cg*cb-db*sg*sb  -dg*sg  0;
         -db*cb               0      0
         ];
    alpha_p = W*[da db dg]'+U*[dda ddb ddg]';

end

function S = skew(v)
    if isvector(v)
        if (length(v)==3)
            % SO(3) case
            S = [  0   -v(3)  v(2)
                  v(3)  0    -v(1)
                 -v(2) v(1)   0];
        elseif (length(v)==2)
            % SO(2) case
            S = [0 -v; v 0];
        else
            error('SMTB:skew:badarg', 'argument must be a 1- or 3-vector');
        end
    else
        error('SMTB:skew:badarg', 'argument must be a 1- or 3-vector');
    end
end
