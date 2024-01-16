clc; clear all;
syms th real;

Sth = sin(th);
Cth = cos(th);

an = [1 -(1-Cth)/th (th-Sth)/th];
bn = [1 (1-Cth)/th (th-Sth)/th];
a = an; a(2:3) = a(2:3)./[th th^2];
b = bn; b(2:3) = b(2:3)./[th th^2];
ain = [1 th/2 1-th*Sth/(2*(1-Cth))];
bin = [1 -th/2 1-th*Sth/(2*(1-Cth))];
ai = ain; ai(2:3) = ai(2:3)./[th th^2];
bi = bin; bi(2:3) = bi(2:3)./[th th^2];
simplify(so3polyprod(an,ain))
simplify(so3polyprod(bn,bin))

%%
da = simplify(diff(a,th)/th);
dda = simplify(diff(da,th)/th);
dai = simplify(diff(ai,th)/th);
ddai = simplify(diff(dai,th)/th);

dan = da; dan(2:3) = dan(2:3).*[th^3 th^4];
ddan = dda; ddan(2:3) = ddan(2:3).*[th^5 th^6];

dain = dai; dain(2:3) = dain(2:3).*[th^3 th^4];
ddain = ddai; ddain(2:3) = ddain(2:3).*[th^5 th^6];

%%
db = simplify(diff(b,th)/th);
ddb = simplify(diff(db,th)/th);
dbi = simplify(diff(bi,th)/th);
ddbi = simplify(diff(dbi,th)/th);

dbn = db; dbn(2:3) = dbn(2:3).*[th^3 th^4];
ddbn = ddb; ddbn(2:3) = ddbn(2:3).*[th^5 th^6];

dbin = dbi; dbin(2:3) = dbin(2:3).*[th^3 th^4];
ddbin = ddbi; ddbin(2:3) = ddbin(2:3).*[th^5 th^6];
