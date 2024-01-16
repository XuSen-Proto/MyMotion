function S = SE3AdsIdx(RphR,i1,i2)
Map = [9 12 15 0 3 6];
iSh = [4 1]; 
iSe = [3 6];
S = zeros(6,1);
S(1:iSe(i1)) = RphR(i2+Map(iSh(i1):6));
end