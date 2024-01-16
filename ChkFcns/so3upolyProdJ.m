function Jo = so3upolyProdJ(c,w,J)
if c(3)~=0
    if c(2)~=0
        c(3) = c(3)/c(2);
    end
    Jo = Ahat3B(c(3)*w,J);
else
    Jo = zeros(size(J));
end
if c(2)~=0
    if c(1)~=0
        c(2) = c(2)/c(1);
    end
    Jo = Ahat3B(c(2)*w,Jo+J);
else
    Jo = Ahat3B(w,Jo);
end
if c(1)~=1
    Jo = c(1)*(Jo+J);
elseif c(1)~=0
    Jo = Jo+J;
end
end