function o = quatxAxis(q,r)
o = [q(2:4) [q([1 4]); -q(3)] [-q(4); q(1:2)] [q(3); -q(2); q(1)]]*r;
end