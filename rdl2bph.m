% Author:     Changmeng Hou(Harbin Engineering University)

% Generate all joints location from platform radius, joint distance and limb length.

function [B,P,H] = rdl2bph(RDL)

assert( ( ( (RDL(3)<=RDL(1)) && (RDL(4)<=RDL(2)) )||(RDL(5)<=RDL(6)) ),...
    'Error Joint Distance!')

rb = RDL(1);
rp = RDL(2);
db = RDL(3);
dp = RDL(4);

etap = pi/3 -asin(dp/2/rp);
etab = asin(db/2/rb);
etaP = [etap 2/3*pi-etap 2/3*pi+etap 4/3*pi-etap 4/3*pi+etap -etap];
etaB = [etab 2/3*pi-etab 2/3*pi+etab 4/3*pi-etab 4/3*pi+etab -etab];
wL = [0 0 0]';
wB = [0 0 0]';
P = [rp*cos(etaP);rp*sin(etaP);[0 0 0 0 0 0]] + wL;
B = [rb*cos(etaB);rb*sin(etaB);[0 0 0 0 0 0]] + wB;


limbMin = RDL(5);
limbMax= RDL(6);

syms h
R_PB = eye(3);
equ1 = norm([0 0 h].'+R_PB*P(:,1)-B(:,1)) == limbMin;
equ2 = norm([0 0 h].'+R_PB*P(:,1)-B(:,1)) == limbMax;
Hmin = double(solve(equ1,h));
Hmax = double(solve(equ2,h));
clear h
if isempty(Hmin), Hmin = 0; end
H = [Hmin;Hmax];

end


