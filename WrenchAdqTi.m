function wo = WrenchAdqTi(wb, tf)
R = quat2R(tf(1:4), 1.0);
f = R*wb(4:6, :);
wo = [R*wb(1:3,:) + hat3(tf(5:7))*f; f];
end