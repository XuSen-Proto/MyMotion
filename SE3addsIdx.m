function S = SE3addsIdx(V,i1,i2)
Map = [2 3; 1 3; 1 2];
im = [1 2 1];
V0 = [V(1:3) V(4:6)];
iSh = [1 2];
map = Map(i2,:);
Stmp = V0(map([2 1]),iSh(i1):-1:1);
Stmp(im(i2),:) = -Stmp(im(i2),:);
S0 = zeros(3,2);
S0(map,1:iSh(i1)) = Stmp;
S = S0(:);
end