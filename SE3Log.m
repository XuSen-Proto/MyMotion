function S = SE3Log(T)
tol = 1e-10;
R = T(1:3,1:3);
p = T(1:3,4);
Cth = (R(1,1)+R(2,2)+R(3,3)-1)/2;
Cth = min(max(Cth,-1),1);
c0 = (1-Cth)<tol;
cp = (Cth+1)<tol;
if (~c0)&&(~cp)
    th = acos(Cth); 
    Sth = sin(th);
    Omg = (th/(2*Sth))*[R(3,2)-R(2,3); R(1,3)-R(3,1); R(2,1)-R(1,2)];
    th_1 = 1/th;
    w2c = (th_1-Sth/(2*(1-Cth)))*th_1;
    v = SE3LogV(Omg,w2c,p);
    S = [Omg; v];
    return;
elseif cp
    th = pi;
    v3 = R(3,3)+1;
    Omg = [R(1:2,3); v3];
    Omg = (th/sqrt(2*v3))*Omg;
    w2c = (1/th)^2;
    v = SE3LogV(Omg,w2c,p);
    S = [Omg; v];
    return;
else
    S = [zeros(3,1);p];
    return;
end
end