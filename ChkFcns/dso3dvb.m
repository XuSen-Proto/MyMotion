function [A,C] = dso3dvb(r,t)
th2 = r'*r;
if th2~=0
    th = sqrt(th2);
    Sth = sin(th);
    Cth = cos(th);
    r = r/th; 
    t = t/th;
    
    % c0 = 1-Cth;
    % c1 = c0/th;
    % c2 = (th-Sth)/th;
    % c34 = [c1 c1^2+c2^2-c2]./(c1^2+(c2-1)^2);
    % cRT = [1 -Sth c0];
    % cA = [1 c34];
    % cCth = so3polyprod(cRT,[-c1 -c2 0]);
    % cCth = so3polyprod(cA,cCth);
    % % cRAT = so3polyprod([1 Sth c0],[1 -c34(1) c34(2)]);
    % c5 = (th*Sth+2*Cth-2)/th;
    % c6 = (3*Sth-2*th-th*Cth)/th;
    % cCtrT = so3polyprod(cRT,[-c5 -c6 0]);
    % cCtrT = so3polyprod(cA,cCtrT);
    % A = so3upoly(cA,r)
    % % RAT = so3upoly(cRAT,r);
    % Cthat = so3upoly(cCth,r)
    % CtrT = so3upoly(cCtrT,r)
    
    c0 = Cth-1;
    c1 = c0/th;
    c2 = th/2;
    c3 = 1+Sth/(2*c1);
    A = so3upoly([1 c2 c3],r);
    Cthat = so3upoly([c1 c3 c1+c2],r);
    CtrT = so3upoly(-[Sth+2*c1 3*c3+0.5*c0 2*c1+Sth/2+c2],r);
    
    Cthat = Ahat3B(t,Cthat');
    rht = Ahat3B(r,t);
    C = (CtrT*rht*r'+Cthat'+Ahat3B((Sth/th-1)*rht,A)')*A;
else
    A = eye(3);
    C = A;
end