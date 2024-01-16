function os = so3TanMap(rs,is,dir,right)
% dir 0: so3 to angular vel; 1: ang_v to so3 
% right 0:wb 1:ws
os = is;
right = (-1)^right;
tmp = zeros(3,1);
if size(is,3) == 1
    for i = 1:size(rs,2)
        th2 = rs(:,i)'*rs(:,i);
        th = sqrt(th2);
        if th<1e-10
            continue;
        end
        sth = sin(th);
        cthm1 = cos(th)-1;
        ru = rs(:,i)/th;
        if dir == 0
            c1 = right*cthm1/th;
            c2 = 1-sth/th;
        else
            c1 = right*0.5;
            th_1 = 1/th;
            cthm1 = 2*cthm1;
            c2 = (th_1 + sth/cthm1)*th_1;
        end
        td = is(:,i,1);
        tmp = c1*td + cross(ru, c2*td);
        os(:,i,1) = td + cross(ru, tmp);
    end
else
    for i = 1:size(rs,2)
        th = sqrt(rs(:,i)'*rs(:,i));
        if th<1e-10
            continue;
        end
        sth = sin(th);
        cthm1 = cos(th)-1;
        th_1 = 1/th;
        ru = rs(:,i)*th_1;
        if dir == 0
            c1 = cthm1*th_1;
            c2 = sth*th_1;
            c3 = -right*(2*c1*th_1 + c2);
            c2 = 1-c2;
            c4 = -(3*c2*th_1 + c1);
            c1 = right*c1;
            id = is(:,i,1);
            rhi = cross(ru, c2*id);
            od = id + cross(ru, c1*id + rhi);
            rd = id;
            rtrd = ru'*rd;
            c3 = rtrd*c3;
            c4 = rtrd*c4;
            rd = rd*th_1;
            od2 = cross(rd, rhi);
            tmp = c3*id + cross(ru, c4*id);
            od2 = od2 + cross(ru, tmp);
        else
            c1 = right*0.5*th;
            cthm1 = 1/(2*cthm1);
            c2 = 1 + sth*cthm1*th;
            c4 = -(2*th_1 + (sth+th)*cthm1);
            id = is(:,i,1);
            rhi = cross(ru, c2*id);
            od = id + cross(ru, c1*id + rhi);
            rd = od;
            rtrd = ru'*rd;
            c4 = rtrd*c4;
            rd = rd.*th_1;
            tmp = cross(rd, id);
            od2 = cross(rd, c1*id + rhi) + c2*cross(ru, tmp);
            tmp = c4*id;
            od2 = od2 + ru*(ru'*tmp)-tmp;
        end
        os(:,i,1) = od;
        id2 = is(:,i,2);
        tmp = c1*id2 + cross(ru, c2*id2);
        od2 = od2 + id2 + cross(ru, tmp);
        os(:,i,2) = od2;
    end
end
end