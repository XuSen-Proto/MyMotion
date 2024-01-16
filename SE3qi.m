function tfi = SE3qi(tf)
tfi = [quati(tf(1:4)); -quat2R(tf(1:4), 1.0)'*tf(5:7)];
end