function q = quat_ZXY(x)
x = x./2;
q = quatx(quatx([cos(x(1));0;0;sin(x(1))],[cos(x(2));sin(x(2));0;0]),[cos(x(3));0;sin(x(3));0]);
end