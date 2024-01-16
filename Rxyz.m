function R = Rxyz(th,i)
Sth = sin(th);
Cth = cos(th);
switch i
    case 1
        it = int8([2 3]);
        r = [1 Cth Cth];
    case 2
        it = int8([3 1]);
        r = [Cth 1 Cth];
    case 3
        it = int8([1 2]);
        r = [Cth Cth 1];
end
% idx = int8([2 3; 3 1; 1 2]);
% it = idx(i,:);
% r = [1 1 1];
% r(it) = Cth;
R = diag(r);
R(it(1),it(2)) = -Sth;
R(it(2),it(1)) = Sth;
end