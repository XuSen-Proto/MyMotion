function [vals,Ans] = My_so3s2wbd(rds)
len = size(rds,2); order = size(rds,3)-1;
vals = zeros(size(rds,1),size(rds,2),order+1);
vals(:,:,1) = rds(:,:,1);
nr = sqrt(sum(rds(:,:,1).*rds(:,:,1),1));
sr = sin(nr); cr = cos(nr);
rds0 = rds; rds = rds./nr;
rtrdn = sum(rds(:,:,1).*rds(:,:,2),1);
if order>2
    rdtrd2n = sum(rds(:,:,1).*rds(:,:,3)+rds(:,:,2).*rds(:,:,2),1);
end
an = cr - 1; bn = nr - sr;
and = -2*cr - nr.*sr + 2;
bnd = 3*sr - nr.*cr - 2*nr;
and2 = nr.*(sr - nr.*cr);
bnd2 = -nr.*and;
if nargout==1
    for i = 1:len
        if nr(:,i)<1e-16
            vals(:,i,2:end) = rds0(:,i,2:end);
        else
            skr = hat3(rds(:,i,1));
            An = eye(3)*nr(i) + (an(i)*eye(3)+bn(i)*skr)*skr;
            vals(:,i,2) = An*rds(:,i,2);
            if order>1
                skrd = hat3(rds(:,i,2)); 
                Aldn = (rtrdn(i)*(and(i)*eye(3)+bnd(i)*skr)*skr + bn(i)*skrd*skr);
                vals(:,i,3) = An*rds(:,i,3)+Aldn*rds(:,i,2);
                if order>2
                    skrd2 = hat3(rds(:,i,3)); ktmp = [rdtrd2n(i) rtrdn(i)^2]';
                    atmp = [and(i) (and2(i)-4*and(i))]*ktmp; btmp = [bnd(i) (bnd2(i)-5*bnd(i))]*ktmp;  
                    vals(:,i,4) = An*rds(:,i,4) + (2*Aldn + an(i)*skrd + bn(i)*skr*skrd)*rds(:,i,3) + ...
                        (atmp*eye(3)+btmp*skr + 2*bnd(i)*rtrdn(i)*skrd + bn(i)*skrd2)*skr*rds(:,i,2);
                end
            end
        end
    end
else
    Ans = zeros(3,3,len);
    for i = 1:len
        if nr(:,i)<1e-16
            vals(:,i,2:end) = rds0(:,i,2:end);
        else
            skr = hat3(rds(:,i,1));
            An = eye(3)*nr(i) + (an(i)*eye(3)+bn(i)*skr)*skr;
            Ans(:,:,i) = An*(1/nr(i));
            vals(:,i,2) = An*rds(:,i,2);
            if order>1
                skrd = hat3(rds(:,i,2)); 
                Aldn = (rtrdn(i)*(and(i)*eye(3)+bnd(i)*skr)*skr + bn(i)*skrd*skr);
                vals(:,i,3) = An*rds(:,i,3)+Aldn*rds(:,i,2);
                if order>2
                    skrd2 = hat3(rds(:,i,3)); ktmp = [rdtrd2n(i) rtrdn(i)^2]';
                    atmp = [and(i) (and2(i)-4*and(i))]*ktmp; btmp = [bnd(i) (bnd2(i)-5*bnd(i))]*ktmp;  
                    vals(:,i,4) = An*rds(:,i,4) + (2*Aldn + an(i)*skrd + bn(i)*skr*skrd)*rds(:,i,3) + ...
                        (atmp*eye(3)+btmp*skr + 2*bnd(i)*rtrdn(i)*skrd + bn(i)*skrd2)*skr*rds(:,i,2);
                end
            end
        end
    end
end
end