function o = quatx21Axis(q,r)
o = [r(2:4) [r(1) -r(4) r(3)]' [r([4 1]); -r(2)] [-r(3); r([2 1])]]*q;
end