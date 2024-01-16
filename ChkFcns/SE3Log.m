function Vhat = SE3Log(T)
% *** CHAPTER 3: RIGID-BODY MOTIONS ***
% Takes a transformation matrix T in SE(3).
% Returns the corresponding se(3) representation of exponential 
% coordinates.
% Example Input:
% 
% clear; clc;
% T = [[1, 0, 0, 0]; [0, 0, -1, 0]; [0, 1, 0, 3]; [0, 0, 0, 1]];
% expmat = SE3Log(T)
% 
% Output:
% expc6 =
%         0         0         0         0
%         0         0   -1.5708    2.3562
%         0    1.5708         0    2.3562
%         0         0         0         0

[R, p] = SE32Rp(T); 
axang = rotm2axang(R);
theta = axang(4);
OmgHat = skew(theta*axang(1:3)');
if isequal(OmgHat, zeros(3))
    Vhat = [zeros(3), T(1: 3, 4); 0, 0, 0, 0];
else
%     theta = acos((trace(R) - 1) / 2);
    Vhat = [OmgHat, (eye(3) - OmgHat / 2 ...
                      + (1 / theta - cot(theta / 2) / 2) ...
                        * OmgHat * OmgHat / theta) * p;
              0, 0, 0, 0];    
end
end