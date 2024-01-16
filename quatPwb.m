function qd = quatPwb(q,wb)
qd = [[-q(2); q([1 4]); -q(3)] [-q(3:4); q(1:2)] [-q(4); q(3); -q(2); q(1)]]*wb;
end