clc; clear all;
syms r t v1 v2 [3 1] real;
syms th0 real;
th = sqrt(r'*r);
rn = r/th; tn = t/th;
v1n = v1/th; v2n = v2/th;
Sth = sin(th);
Cth = cos(th);
rnh = hat3(rn);
R = eye(3)+Sth*rnh+(1-Cth)*rnh*rnh;
A = eye(3)-(1-Cth)/th*rnh+(th-Sth)/th*rnh*rnh;
%%
Ai = eye(3)+th/2*rnh+((th*sin(th))/(2*cos(th)-2)+1)*rnh*rnh;
an = [1 th0/2 1-th0*sin(th0)/(2*(1-cos(th0)))];
a = an; a(2:3) = a(2:3)./[th0 th0^2];
da0 = simplify(diff(a,th0)/th0);
dda0 = simplify(diff(da0,th0)/th0);
dda = subs(dda0,th);
da = subs(da0,th);
dan = da; dan(2:3) = dan(2:3).*[th^3 th^4];
ddan = dda; ddan(2:3) = ddan(2:3).*[th^5 th^6];
Aiv1 = Ai*v1;
dAiv1dr = simplify(jacobian(Aiv1,r));
dAv1drTv2 = dAiv1dr'*v2;
dAv1drTv2_dr = simplify(jacobian(dAv1drTv2,r));
dAv1drTv2_dv1 = simplify(jacobian(dAv1drTv2,v1));
matlabFunction(dAv1drTv2_dr,dAv1drTv2_dv1,'File','dso3dwbchk','Vars',{r,v1,v2});
%%
matlabFunction(dAiv1dr,dan,ddan,'File','dAiv1drchk','Vars',{r,v1});
%%
an = [1 -(1-cos(th0))/th0 (th0-sin(th0))/th0];
a = an; a(2:3) = a(2:3)./[th0 th0^2];
da0 = simplify(diff(a,th0)/th0);
dda0 = simplify(diff(da0,th0)/th0);
dda = subs(dda0,th);
da = subs(da0,th);
C = so3upoly(da,r);
% v2TCv1 = v1'*C'*v2;
% Cm = simplify(jacobian(C'*v2,r));

%%
dan = da; dan(2:3) = dan(2:3).*[th^3 th^4];
ddan = dda; ddan(2:3) = ddan(2:3).*[th^5 th^6];
rhv = Ahat3B(rn,v1n);
dA = hat3(dan(3)*rn);
dA(1:4:9) = dA(1:4:9)+dan(2);
C = dA;
dArhv = dA*rhv;
dA = dArhv*rn';
Cthat = hat3(-an(3)*rn);
Cthat(1:4:9) = Cthat(1:4:9)-an(2);
Cthat = Cthat*hat3(v1n);
dA = dA + Cthat - dan(3)*hat3(rhv);
% Am = simplify(rn*(v2/th)'*dA);
%%
Cm = simplify(r*v1'*jacobian(C'*v2,r));
Am = simplify(rn*(v1/th)'*dA);
matlabFunction(v2TCv1,Cm*th^2,C*th^2,dan,ddan,'File','Cmchk','Vars',{r,v1,v2});
%%
ATv1 = A'*v1;
dATv1_dr = simplify(jacobian(ATv1,r));
matlabFunction(dATv1_dr,'File','dATv1_drchk','Vars',{r,v1});
%%
Av1 = A*v1;
dAv1_dr = simplify(jacobian(Av1,r));
dAv1drTv2 = dAv1_dr*v2;
dAv1drTv2_dr = simplify(jacobian(dAv1drTv2,r));
dAv1drTv2_dv1 = simplify(jacobian(dAv1drTv2,v1));
an1 = subs(an,th);
dAdv = hat3(an1(2)*v2n)+an1(3)*(hat3(rn)*hat3(v2n)+hat3(v2n)*hat3(rn))+rn'*v2n*C*hat3(rn);
%%
matlabFunction(dAv1drTv2_dr,dAv1drTv2_dv1,C,'File','dwbdso3chk','Vars',{r,v1,v2});
%%
p = t+(1-Cth)*rnh*tn+(th-Sth)*rnh*rnh*tn;
dp = jacobian(p,[r;t]);
dvb = simplify(R'*dp);
dvbf = matlabFunction(dvb,'Vars',{r,t});
Af = matlabFunction(A,'Vars',{r});
%%
clc
rt = rand(3,1);
tt = rand(3,1);
%%
[A1, C1] = dvbdso3(rt,tt);
[A2, C2] = dso3dvb(rt,tt);
[A1 zeros(3); C1 A1]*[A2 zeros(3); C2 A2]
% dvbf(rt,tt)-[C1 A1]

%%
clc; clear all;
syms th real;
Sth = sin(th);
Cth = cos(th);
c0 = 1-Cth;
c2 = (th-Sth)/th;
cRT = [1 -Sth c0];
c1 = c0/th;

cCth = simplify(so3polyprod(cRT,[c1 c2 0]));
c3 = (th*Sth+2*Cth-2)/th;
c4 = (3*Sth-2*th-th*Cth)/th;
cCtrT = simplify(so3polyprod(cRT,[c3 c4 0]));
pretty(cCth)
pretty(cCtrT)

%%
clc; clear all;
syms th real;
Sth = sin(th);
Cth = cos(th);
c0 = 1-Cth;
c1 = c0/th;
c2 = (th-Sth)/th;
c34 = simplify([c1 c1^2+c2^2-c2]./(c1^2+(c2-1)^2));
cRT = [1 -Sth c0];
cA = [1 c34];
cCth = so3polyprod(cRT,[-c1 -c2 0]);
cCth = simplify(so3polyprod(cA,cCth));
cRAT = simplify(so3polyprod([1 Sth c0],[1 -c34(1) c34(2)]));
c5 = -(th*Sth+2*Cth-2)/th;
c6 = -(3*Sth-2*th-th*Cth)/th;
cCtrT = so3polyprod(cRT,[c5 c6 0]);
cCtrT = simplify(so3polyprod(cA,cCtrT));
pretty(c34)
pretty(cCth)
pretty(cCtrT)
pretty(cRAT)
%%
c = hat3(r)*t;
matlabFunction(c,'File','Ahat3B','Vars',{r,t});