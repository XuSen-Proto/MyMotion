function wo = WrenchAdqT(wb, tf)
    R = quat2R(tf(1:4), 1.0)';
    wo = [R*(wb(1:3) - cross(tf(5:7),wb(4:6))); R*wb(4:6)];
end