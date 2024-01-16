clc; clear all;
r = rand(3,1)-0.5;
t = rand(3,1)-0.5;
th = norm(r);
rn = r/th;
R1 = SO3Exp(rn,th);
[A,C] = dvbdso3(r,t);
[Ai,Ci] = dso3dvb(r,t);
Ro = SO3Exps(R1,th);
A1 = Ro.dwbdso3(eye(3));
C1 = Ro.d_dwbdso3Pt(t);
A1i = Ro.dso3dwb(eye(3));
C1i = Ro.d_dso3dwbPt(t);
max(abs([A-A1 C-C1 A1i-Ai C1i-Ci]),[],'all');

%%
clc; clear all;
r = rand(3,1)-0.5;
t = rand(3,1)-0.5;
v = rand(3,1)-0.5;
th = norm(r);
tn = t/th;
vn = v/th;
T = SE3Exp([r;t],[]);
To = SE3s([r;t],[]);
To.GetRTp();
[c1,c2,c3] = To.dwdrCoefs(1);
To.getA(c1);
A = To.A;
dA = To.dAvn_dr(tn,c1(2:3),c2);
% [c1,c2,c3] = To.drdwCoefs(1);
% To.getA(c1);
% Ai = To.A;
% dAi = To.dAvn_dr(tn,c1(2:3),c2);
% [dAic,dan,ddan] = dAiv1drchk(r,t);

% A1 = dATv1_drchk(r,t);
% A2 = To.dAvn_dr(t/th,[-To.c1(2) To.c1(3)],[-To.c2(1) To.c2(2)]);

% [val,Cm,C,da,dda] = Cmchk(r,t,v)
[dA1c,dA2c,C] = dwbdso3chk(r,t,v)
[dA1,dA2] = To.ddAvndrmn(tn,vn,c1(2:3),c2,c3);


%%
r = rand(3,1)-0.5;
t = rand(3,1)-0.5;
th = norm(r);
rn = r/th;
Sth = sin(th);
Cth = cos(th);
c0 = 1-Cth;
f1 = @()Rt1_mex(rn,t,Sth,c0);
f2 = @()Rt2_mex(rn,t,Sth,c0);
timeit(f1)
timeit(f2)
%%
for i = 1:20000
t1 = Rt1_mex(rn,t,Sth,c0);
t2 = Rt2_mex(rn,t,Sth,c0);
end

%%
cfg = coder.config('mex');
cfg.TargetLang = 'C++';
cfg.EchoExpressions = false;
cfg.SaturateOnIntegerOverflow = false;
cfg.ResponsivenessChecks = false;
codegen -config cfg Rt2 -args {zeros(3,1),zeros(3,1),0,0} -report -O enable:inline -O enable:openmp

%%
function t = Rt1(rn,t,Sth,c0)
t = t+Ahat3B(rn,Sth*t+Ahat3B(rn,c0*t));
end

%%
function t = Rt2(rn,t,Sth,c0)
t = so3upoly([1 Sth c0],rn)*t;
end