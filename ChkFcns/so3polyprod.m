function c = so3polyprod(c1,c2)
C = c1'*c2;
c = [C(1) C(2)+C(4)-C(6)-C(8) C(3)+C(5)+C(7)-C(9)];
end