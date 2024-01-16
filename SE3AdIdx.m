function S = SE3AdIdx(RphR,i1,i2)
Map = int8([-2 -1 0 7 8 9]);
iSh = int8([1 4]); 
iSe = int8([6 3]);
S = zeros(6,1);
S(iSh(i1):6) = RphR(3*i2+Map(1:iSe(i1)));
end