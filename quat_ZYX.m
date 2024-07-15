function q = quat_ZYX(x)
x = x./2;
q = quatx(quatx([cos(x(1));0;0;sin(x(1))],[cos(x(2));0;sin(x(2));0]),[cos(x(3));sin(x(3));0;0]);
end