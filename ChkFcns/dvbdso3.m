function [A,C] = dvbdso3(r,t)
th2 = r'*r;
th = sqrt(th2);
Sth = sin(th);
Cth = cos(th);
r = r/th; 
t = t/th;
c0 = 1-Cth;
c1 = c0/th;
c2 = (th-Sth)/th;
c3 = 2*c1-Sth;

% cRT = [1 -Sth c0];
% cCth = so3polyprod(cRT,[c1 c2 0]);
% c4 = (th*Sth+2*Cth-2)/th;
% c5 = (3*Sth-2*th-th*Cth)/th;
% cCtrT = so3polyprod(cRT,[c4 c5 0]);
% Cthat = Ahat3B(t,Cthat');
% CtrT = so3upoly(cCtrT,r);

rht = Ahat3B(r,t);
R = so3upoly([1 Sth c0],r);
A = so3upoly([1 -c1 c2],r);
Cthat = so3upoly([c1 c2-c0 c3],r);
Cthat = Ahat3B(t,Cthat');
CtrT = so3upoly([-c3 Sth*c1-(2*Cth+1)*c2 Sth*(2+c2)-4*c1],r);
C = CtrT*rht*r'+Cthat'+Ahat3B(c2*rht,R)';
% J = [A zeros(3); C A];
end