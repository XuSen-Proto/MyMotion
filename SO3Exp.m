function R = SO3Exp(Omg,thetain)
if isempty(thetain)
    theta = norm(Omg);
    if abs(theta)<1e-12
        R = eye(3); return;
    end
%     Omg = Omg/theta;
%     R = SO3Rodrigues(Omg,theta);
    R = so3upoly([1 sin(theta)/theta (1-cos(theta))/(theta^2)],Omg);
    return;
else
    theta = thetain;
%     R = SO3Rodrigues(Omg,theta);
    R = so3upoly([1 sin(theta) 1-cos(theta)],Omg);
    return;
end
end