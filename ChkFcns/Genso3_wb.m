clc; clear all;
syms r wb wbd rd rd2 [3 1] real;
nr = sqrt(r'*r);
nr_1 = 1/nr;
skr = skew(r);
rd2wb = eye(3) - nr_1^2*skr*((1-cos(nr))*eye(3) - (1-nr_1*sin(nr))*skr);
wb2rd = inv(rd2wb);

%% rd 2 wb
wbn = rd2wb*rd;
rd2wb_l = rd2wb(:);
Jrd2wb_l = jacobian(rd2wb_l,r);
wbdn = rd2wb*rd2 + reshape(Jrd2wb_l*rd,3,3)*rd;

%% wb 2 rb
rdn = wb2rd*wb;
wb2rd_l = wb2rd(:);
Jwb2rd_l = jacobian(wb2rd_l,r);
rd2n = wb2rd*wbd + reshape(Jwb2rd_l*rdn,3,3)*wb;

%%
matlabFunction(reshape([wbn wbdn],3,1,[]),'File','so3c2wb','Vars',{r,rd,rd2});
matlabFunction(reshape([rdn rd2n],3,1,[]),'File','wb2so3c','Vars',{r,wb,wbd});
%%
rws = rand(3,1,3);
tic
vals1 = My_wbd2so3(rws);
t1 = toc;
tic
vals2 = wb2so3c(rws(:,:,1),rws(:,:,2),rws(:,:,3));
t2 = toc;
[t1 t2 max(abs(vals1(:,:,2:3) - vals2),[],'all')]
