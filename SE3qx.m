function tfo = SE3qx(tf1, tf2)
tfo = [quatx(tf1(1:4),tf2(1:4)); tf1(5:7)+quat2R(tf1(1:4), 1.0)*tf2(5:7)];
end