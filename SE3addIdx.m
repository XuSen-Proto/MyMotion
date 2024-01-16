function S = SE3addIdx(V,i1,i2)
Map = int8([2 3; 1 3; 1 2]);
im = int8([2 1 2]);
V0 = [V(1:3) V(4:6)];
iSh = 3-i1;
map = Map(i2,:);
im2 = im(i2);
Stmp = V0(map([2 1]),1:iSh);
Stmp(im2,:) = -Stmp(im2,:);
S0 = zeros(3,2);
S0(map,i1:2) = Stmp;
S = S0(:);
end