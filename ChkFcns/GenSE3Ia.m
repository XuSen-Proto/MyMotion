clc; clear all;
syms I0_ [10 1] real;
syms SymMtx6Vec [21 1] real;
syms T1_ [4 4] real;
syms v1_ [6 1] real;
syms D real;
syms R phR [3 3] real;
Ic = zeros(6,6,"like",I0_1);
SymMtx6 = zeros(6,6,"like",I0_1);
Ic(1,2:3) = I0_(4:5);
Ic(2,3) = I0_(6);
Ic(1:3,4:6) = skew(I0_(8:10));
Ic = Ic+Ic'+diag([I0_(1:3); I0_(7)*ones(3,1)]);
Ia_c = [diag(Ic); Ic(1,2:3)'; Ic(2,3); Ic(4,5:6)'; Ic(5,6); reshape(Ic(1:3,4:6),[],1)];
SymMtx6(1,2:3) = SymMtx6Vec(7:8)';
SymMtx6(2,3) = SymMtx6Vec(9);
SymMtx6(4,5:6) = SymMtx6Vec(10:11)';
SymMtx6(5,6) = SymMtx6Vec(12);
SymMtx6(1:3,4:6) = reshape(SymMtx6Vec(13:21),3,[]);
SymMtx6 = SymMtx6+SymMtx6'+diag(SymMtx6Vec(1:6));
R1_ = T1_(1:3,1:3); p1_ = T1_(1:3,4);
Adov_T = [R1_ zeros(3); skew(p1_)*R1_ R1_];
Iacv = SymMtx6*v1_;
IacT = Adov_T'*SymMtx6*Adov_T;
IacTVec = [diag(IacT); IacT(1,2:3)'; IacT(2,3); IacT(4,5:6)'; IacT(5,6); reshape(IacT(1:3,4:6),[],1)];
Iaa = SymMtx6+D*(v1_*v1_');
IaaVec = [diag(Iaa); Iaa(1,2:3)'; Iaa(2,3); Iaa(4,5:6)'; Iaa(5,6); reshape(Iaa(1:3,4:6),[],1)];

AdT = [R zeros(3); phR R];
IcT2 = AdT'*SymMtx6*AdT;
IacTVec2 = [diag(IcT2); IcT2(1,2:3)'; IcT2(2,3); IcT2(4,5:6)'; IcT2(5,6); reshape(IcT2(1:3,4:6),[],1)];
%%
prefix = '../SE3'; 
matlabFunction(Iacv,'File',[prefix 'Iavd'],'Vars',{SymMtx6Vec,v1_});
matlabFunction(IacTVec,'File',[prefix 'IacT'],'Vars',{T1_,SymMtx6Vec});
matlabFunction(Ia_c,'File',[prefix 'Ic2Ia'],'Vars',{I0_});
matlabFunction(IaaVec,'File',[prefix 'IaPlusaaT'],'Vars',{SymMtx6Vec,v1_,D});
matlabFunction(IacTVec2,'File',[prefix 'IacT2'],'Vars',{R,phR,SymMtx6Vec});
%  prefix2 = 'ov_';
% matlabFunction(Iacv,'File',[prefix prefix2 'Iavd'],'Vars',{SymMtx6Vec,v1_});
% matlabFunction(IacTVec,'File',[prefix prefix2 'IacT'],'Vars',{T1_,SymMtx6Vec});
% matlabFunction(Ia_c,'File',[prefix prefix2 'Ic2Ia'],'Vars',{I0_});
% matlabFunction(IaaVec,'File',[prefix prefix2 'IaPlusaaT'],'Vars',{SymMtx6Vec,v1_,D});
% matlabFunction(IacTVec2,'File',[prefix prefix2 'IacT2'],'Vars',{R,phR,SymMtx6Vec});
%%
T = [eul2rotm(rand(1,3)) rand(3,1); 0 0 0 1];
D = 5;
AdT = SE3Adm(T);
Ic = rand(10,1);
vd = rand(6,1);
Ia = SE3Ic2Ia(Ic);
Ia2v = rand(21,1);
Ia2 = zeros(6,6);
Ia2(1,2:3) = Ia2v(7:8)';
Ia2(2,3) = Ia2v(9);
Ia2(4,5:6) = Ia2v(10:11)';
Ia2(5,6) = Ia2v(12);
Ia2(1:3,4:6) = reshape(Ia2v(13:21),3,[]);
Ia2 = Ia2+Ia2'+diag(Ia2v(1:6));
Icvd = SE3Ivd(Ic,vd);
Iavd = SE3Iavd(Ia,vd);
Ia2vd1 = SE3Iavd(Ia2v,vd);
Ia2vd2 = Ia2*vd;
Ia3v = SE3IacT(T,Ia2v);
Ia3 = AdT'*Ia2*AdT;
Ia3vr = [diag(Ia3); Ia3(1,2:3)'; Ia3(2,3); Ia3(4,5:6)'; Ia3(5,6); reshape(Ia3(1:3,4:6),[],1)];
Ia4v = SE3IaPlusaaT(Ia2v,vd,D);
Ia4 = Ia2+(D*vd)*vd';
Ia4vr = [diag(Ia4); Ia4(1,2:3)'; Ia4(2,3); Ia4(4,5:6)'; Ia4(5,6); reshape(Ia4(1:3,4:6),[],1)];

[max(abs(Icvd-Iavd)) max(abs(Ia2vd1-Ia2vd2)) max(abs(Ia3v-Ia3vr)) max(abs(Ia4v-Ia4vr))]
