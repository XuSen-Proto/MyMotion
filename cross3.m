function [o,do1,do2] = cross3(i1,i2,di1,di2)
% R3的叉乘及其微分映射
% 第三维度表示阶数，两阶以上第二维度只能是1
o = cross(i1,i2);
if nargout==2
    if nargin == 3
        do1 = -pagemtimes(hat3(i2),di1);
    elseif nargin == 4
        % do1 = pagemtimes(hat3(i1),di2);
        for i = 1:size(di1,3)
            do1(:,:,i) = hat3(i1)*di2(:,:,i);
        end
        if ~isempty(di1)
            % do1 = do1 - pagemtimes(hat3(i2),di1);
            for i = 1:size(di1,3)
                do1(:,:,i) = do1(:,:,i) - hat3(i2)*di1(:,:,i);
            end
            if size(di1,3)==2
                do1(:,:,2) = do1(:,:,2) + 2*cross(di1(:,:,1),di2(:,:,1));
            end
        end
    end
elseif nargout == 3
    do1 = -pagemtimes(hat3(i2),di1);
    do2 = pagemtimes(hat3(i1),di2);
    if size(di1,3)==2
        do1(:,:,2) = do1(:,:,2) + cross(di1(:,:,1),di2(:,:,1));
        do2(:,:,2) = do2(:,:,2) + cross(di1(:,:,1),di2(:,:,1));
    end
end
end