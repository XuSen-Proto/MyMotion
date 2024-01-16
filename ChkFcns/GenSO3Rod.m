clc;clear all;
syms Omega v pn [3 1] real;
syms theta w2c real;
OmgHat = skew(Omega);
Sth = sin(theta); Vth = 1-cos(theta);
R = eye(3)+(Sth*eye(3)+Vth*OmgHat)*OmgHat;
ov = cross(Omega,v); oov = cross(Omega,ov);
p = theta*v+Vth*ov+(theta-Sth)*oov;
T = [R p; 0 0 0 1];
% ctg_th2 = (1+c_theta)/s_theta;
vp = (eye(3)-OmgHat/2+w2c*OmgHat*OmgHat)*pn;
%%
matlabFunction(R,'File','SO3Rodrigues','Vars',{Omega,theta});
matlabFunction(T,'File','SE3Rodrigues','Vars',{Omega,theta,v});
matlabFunction(vp,'File','SE3LogV','Vars',{Omega,w2c,pn});
%%
Omega = rand(3,1);
theta = norm(Omega);
Omega = Omega/theta;    
v = rand(3,1);
tic
OmgHat = skew(Omega);
Sth = sin(theta); Vth = 1-cos(theta);
R1 = eye(3)+(Sth*eye(3)+Vth*OmgHat)*OmgHat;
t1 = toc;
tic
R2 = SO3Rodrigues(Omega,theta);
t2 = toc;

tic
OmgHat = skew(Omega);
Sth = sin(theta); Vth = 1-cos(theta);
R = eye(3)+(Sth*eye(3)+Vth*OmgHat)*OmgHat;
ov = cross(Omega,v); oov = cross(Omega,ov);
p = theta*v+Vth*ov+(theta-Sth)*oov;
T1 = [R p; 0 0 0 1];
t3 = toc;
tic
T2 = SE3Rodrigues(Omega,theta,v);
t4 = toc;

thetan = pi;
Rn = SO3Rodrigues(Omega,thetan);
tic
axang = rotm2axang(Rn);
t5 = toc;
tic
[Omg,th,Sth,Cth] = SO3Log(Rn);
t6 = toc;

pn = rand(3,1);
T = [R pn; 0 0 0 1];

tic
S1h = SE3Log(T);
t7 = toc;
S1 = [vex(S1h(1:3,1:3)); S1h(1:3,4)];

tic
S2 = SE3Log(T);
t8 = toc;

[max(abs(R1-R2),[],'all') max(abs(T1-T2),[],'all')]
[t1 t2]
[t3 t4]
[t5 t6]
[t7 t8]

