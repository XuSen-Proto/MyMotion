function M = so3upoly(c,w)
M = hat3(c(3)*w);
M(1:4:9) = M(1:4:9)+c(2);
M = Ahat3B(w,M);
M(1:4:9) = M(1:4:9)+c(1);
end