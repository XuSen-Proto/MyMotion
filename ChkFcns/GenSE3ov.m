clc; clear all;
syms R phR [3 3] real;
syms RphR [3 6] real;
% syms p1_ p2_ [3 1] real;
syms T1_ T2_ [4 4] real;
syms v1_ v2_ f [6 1] real;
syms I0_ [10 1] real;
% syms Ic(I0_) [6 6] real;
Ic = zeros(6,6,"like",v1_1);
Ic(1,2:3) = I0_(4:5);
Ic(2,3) = I0_(6);
Ic(1:3,4:6) = skew(I0_(8:10));
Ic = Ic+Ic'+diag([I0_(1:3); I0_(7)*ones(3,1)]);
% T1 = [R1_ p1_; 0 0 0 1];
% T2 = [R2_ p2_; 0 0 0 1];

R1_ = T1_(1:3,1:3); p1_ = T1_(1:3,4);
R2_ = T2_(1:3,1:3); p2_ = T2_(1:3,4);

T1T2 = [R1_*R2_ p1_+R1_*p2_; 0 0 0 1];
Adov_T = [R1_ zeros(3); skew(p1_)*R1_ R1_];
adov_v = [skew(v1_(1:3)) zeros(3); skew(v1_(4:6)) skew(v1_(1:3))];

Ad_ov = Adov_T*v1_;
Ads_ov = Adov_T'*f;
ad_ov = adov_v*v2_;
ads_ov = adov_v'*f;
Icv = Ic*v1_;
IcT = Adov_T'*Ic*Adov_T;
IcTc = [diag(IcT(1:3,1:3)); IcT(1,2:3)'; IcT(2,3); IcT(4,4); IcT(3,5); IcT(1,6); IcT(2,4)];

I_1_v = Ic\f;

AdT = [R zeros(3); phR R];
AdTn = [RphR(:,1:3) zeros(3); RphR(:,4:6) RphR(:,1:3)];
AdTns = AdTn';
IcT2 = AdT'*Ic*AdT;
IcTc2 = [diag(IcT2(1:3,1:3)); IcT2(1,2:3)'; IcT2(2,3); IcT2(4,4); IcT2(3,5); IcT2(1,6); IcT2(2,4)];
%%
prefix = '../SE3'; 
matlabFunction(T1T2,'File',[prefix 'x'],'Vars',{T1_,T2_});
matlabFunction(Ad_ov,'File',[prefix 'Ad'],'Vars',{T1_,v1_});
matlabFunction(Ads_ov,'File',[prefix 'Ads'],'Vars',{T1_,f});
matlabFunction(ad_ov,'File',[prefix 'add'],'Vars',{v1_,v2_});
matlabFunction(ads_ov,'File',[prefix 'adds'],'Vars',{v1_,f});
matlabFunction(Icv,'File',[prefix 'Ivd'],'Vars',{I0_,v1_});
matlabFunction(IcTc,'File',[prefix 'IcT'],'Vars',{T1_,I0_});
matlabFunction(I_1_v,'File',[prefix 'IcLdf'],'Vars',{I0_,f});
matlabFunction(IcTc2,'File',[prefix 'IcT2'],'Vars',{R,phR,I0_});
% prefix2 = 'ov_';
% matlabFunction(T1T2,'File',[prefix 'x'],'Vars',{T1_,T2_});
% matlabFunction(Ad_ov,'File',[prefix prefix2 'Ad'],'Vars',{T1_,v1_});
% matlabFunction(Ads_ov,'File',[prefix prefix2 'Ads'],'Vars',{T1_,f});
% matlabFunction(ad_ov,'File',[prefix prefix2 'add'],'Vars',{v1_,v2_});
% matlabFunction(ads_ov,'File',[prefix prefix2 'adds'],'Vars',{v1_,f});
% matlabFunction(Icv,'File',[prefix prefix2 'Ivd'],'Vars',{I0_,v1_});
% matlabFunction(IcTc,'File',[prefix prefix2 'IcT'],'Vars',{T1_,I0_});
% matlabFunction(I_1_v,'File',[prefix prefix2 'IcLdf'],'Vars',{I0_,f});
% matlabFunction(IcTc2,'File',[prefix prefix2 'IcT2'],'Vars',{R,phR,I0_});



