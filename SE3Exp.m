function T = SE3Exp(Scr,thetain)
Omg = Scr(1:3);
v = Scr(4:6); 
if isempty(thetain)
    theta = norm(Omg);
    if theta<1e-10
        T = [eye(3) v; 0 0 0 1]; return;
    end
    Sth = sin(theta);
    Cth = cos(theta);
    Omg = Omg/theta;
%     v = v/theta;
%     T = SE3Rodrigues(Omg,theta,v);
    c0 = 1-Cth;
    c1 = c0/theta;
    c2 = (theta-Sth)/theta;
    R = so3upoly([1 Sth c0],Omg);
    B = so3upoly([1 c1 c2],Omg);
    T = [R B*v; 0 0 0 1];
    return;
else
    theta = thetain;
%     T = SE3Rodrigues(Omg,theta,v);
    Sth = sin(theta);
    Cth = cos(theta);
    c0 = 1-Cth;
    c2 = theta-Sth;
    R = so3upoly([1 Sth c0],Omg);
    B = so3upoly([theta c0 c2],Omg);
    T = [R B*v; 0 0 0 1];
    return;
end
end