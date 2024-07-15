function [vo,dvo] = UnitVec(v,dv)
% 单位化映射及其切映射
% 第三维度表示阶数，两阶以上第二维度只能是1
if nargin == 1
    vo = v./(sqrt(sum(v.*v)));
else
    if nargout == 2
        nv = sqrt(v'*v);
        vo = v./nv;
        dvo = dv./nv;
        
        if size(dv,3) == 1
            dvo = dvo - vo*(vo'*dvo);
        else
            dvo(:,:,2) = dvo(:,:,2) - vo*(vo'*dvo(:,:,2));
            vTdv = vo'*dvo(:,:,1);
            dv2 = dvo(:,:,1)'*dvo(:,:,1);
            dvn = dvo(:,:,1) - vo*(vo'*dvo(:,:,1));
            dvo(:,:,2) = dvo(:,:,2) - vTdv*dvn - dvo(:,:,1)*vTdv + (2*vTdv^2-dv2)*vo;
            dvo(:,:,1) = dvn;
        end
    else
        vo = v./(sqrt(sum(v.*v)));
    end
end