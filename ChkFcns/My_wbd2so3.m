function vals = My_wbd2so3(rws)
len = size(rws,2); order = size(rws,3)-1;
vals = zeros(size(rws,1),size(rws,2),order+1);
vals(:,:,1) = rws(:,:,1);
nr = sqrt(sum(rws(:,:,1).*rws(:,:,1),1));
sr = sin(nr); cr = cos(nr);
% rn = rws(:,:,1)./(ones(3,1)*nr);
an = cr - 1; bn = nr - sr;
and = -2*cr - nr.*sr + 2;
bnd = 3*sr - nr.*cr - 2*nr;
and2 = nr.*(sr - nr.*cr);
% bnd2 = nr.*(2*cr + nr.*sr - 2);
bnd2 = -nr.*and;
for i = 1:len
    if nr(:,i)<1e-16
        vals(:,i,2:end) = rws(:,i,2:end);
    else
        rn = rws(:,i,1)/nr(i);
        skr = hat3(rn);
        An = eye(3)*nr(i) + (an(i)*eye(3)+bn(i)*skr)*skr;
        rdn = An\rws(:,i,2); 
        rtrdn = rn'*rdn; 
        skrd = hat3(rdn);
        vals(:,i,2) = rdn*nr(i);
        if order>1
            Aldn = (rtrdn*(and(i)*eye(3)+bnd(i)*skr)*skr + bn(i)*skrd*skr);
            rd2n = An\(rws(:,i,3) - Aldn*rdn); skrd2 = hat3(rd2n);
            vals(:,i,3) = rd2n*nr(i);
            if order>2
                rdtrd2n = rn'*rd2n+rdn'*rdn; 
                atmp = and(i)*rdtrd2n+(and2(i)-4*and(i))*rtrdn^2;
                btmp = bnd(i)*rdtrd2n+(bnd2(i)-5*bnd(i))*rtrdn^2;
                rd3n = An\(rws(:,i,4) - (2*Aldn + an(i)*skrd + bn(i)*skr*skrd)*rd2n ...
                    -(atmp*eye(3)+btmp*skr + 2*bnd(i)*rtrdn*skrd + bn(i)*skrd2)*skr*rdn);     
                vals(:,i,4) = rd3n*nr(i);
            end
        end
    end
end
end