function W = quat_dRv(h,nh2,v)
h = (2/nh2)*h;
% W = hat3(v)*[h(2:4) [-h(1); h(4); -h(3)] [-h([4 1]); h(2)] [h(3); -h([2 1])]];
W = Ahat3B(v,[h(2:4) [-h(1); h(4); -h(3)] [-h([4 1]); h(2)] [h(3); -h([2 1])]]);
end