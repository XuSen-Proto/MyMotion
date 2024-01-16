clc; clear all;
I0_ = rand(10,1);
Ia_ = rand(21,1);
R = eul2rotm(rand(1,3));
p = rand(3,1);
T = [R p; 0 0 0 1];
v = rand(6,1);
AdT0 = SE3Adm(T);
adv0 = [skew(v(1:3)) zeros(3); skew(v(4:6)) skew(v(1:3))];
phR = skew(p)*R;
Ic0 = zeros(6,6);
Ic0(1,2:3) = I0_(4:5);
Ic0(2,3) = I0_(6);
Ic0(1:3,4:6) = skew(I0_(8:10));
Ic0 = Ic0+Ic0'+diag([I0_(1:3); I0_(7)*ones(3,1)]);
Ia0 = zeros(6,6);
Ia0(1,2:3) = Ia_(7:8)';
Ia0(2,3) = Ia_(9);
Ia0(4,5:6) = Ia_(10:11)';
Ia0(5,6) = Ia_(12);
Ia0(1:3,4:6) = reshape(Ia_(13:21),3,[]);
Ia0 = Ia0+Ia0'+diag(Ia_(1:6));

Ic = zeros(6); Ia = Ic; AdT = Ic; AdTs = Ic; adv=Ic; advs = Ic;
for iS = 1:6
    i1 = ceil(iS/3); i2 = iS-3*(i1-1);
    Ic(:,iS) = SE3IcIdx(iS,I0_);
    Ia(:,iS) = SE3IaIdx(iS,Ia_);
    AdT(:,iS) = SE3AdIdx([R,phR],i1,i2);
    AdTs(:,iS) = SE3AdsIdx([R,phR],i1,i2);
    adv(:,iS) = SE3addIdx(v,i1,i2);
    advs(:,iS) = SE3addsIdx(v,i1,i2);
end
[max(abs(Ic-Ic0),[],'all') max(abs(Ia-Ia0),[],'all') max(abs(AdT-AdT0),[],'all')...
    max(abs(AdTs-AdT0'),[],'all') max(abs(adv-adv0),[],'all') max(abs(advs-adv0'),[],'all')]
