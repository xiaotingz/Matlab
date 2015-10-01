a = eye(3);
b = [0,-1,1;1,0,-1;-1,1,0];
c = b^2;

% NOTICE, this is for P=1, if P=-1, RotationMatrix is its transpose (which hemisphere)
RotatMatrix = a + 1/2*b + 1/6*c;

para = 1/sqrt(2/3);
phi_1 = atan2(para*RotatMatrix(3,1),para*RotatMatrix(3,2));
PHI = acos(RotatMatrix(3,3));
phi_2 = atan2(para*RotatMatrix(1,3),para*RotatMatrix(2,3));

% from one rotation matrix, can always get two sets of euler angle. plus we don't know about P, list the full answer shit
AnsSet = [phi_1,cos(PHI),phi_2;pi-phi_1,cos(pi-PHI),pi-phi_2;0,0,0;phi_2,cos(PHI),phi_1;pi-phi_2,cos(pi-PHI),pi-phi_1];

% the domain in consideration is one eighth sphere
Domain = [pi/2,1,pi/2;0,0,0;0,0,0];
