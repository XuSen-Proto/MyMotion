function [Omg,th,Sth,Cth] = SO3Log(R)
tol = 1e-10;
Cth = (R(1,1)+R(2,2)+R(3,3)-1)/2;
Cth = min(max(Cth,-1),1);
if (1-Cth)<tol
    Omg = zeros(3,1);
    th = 0;
    Sth = 0;
    Cth = 1;
    return;
elseif (Cth+1)<tol
    th = pi;
    Sth = 0;
    Cth = -1;
    Omg = [R(1:2,3); R(3,3)+1];
    Omg = (1/sqrt(2*Omg(3)))*Omg;
else
    th = acos(Cth); 
    Sth = sin(th);
    Omg = (1/(2*Sth))*[R(3,2)-R(2,3); R(1,3)-R(3,1); R(2,1)-R(1,2)];
end
if nargout==1
    Omg = Omg*th;
end
end