clc; clear all;
q = rand(1,4); v = rand(1,3);
tic
R = quat2rotm(q);
v2 = R*v';
toc
tic
qinv = quatinv(q);
v1 = quatmultiply(q,quatmultiply([0 v],qinv));
toc
%%
clc; clear all;
syms q qd [1 4] real;
R = eye(3)+2/(q*q')*[-q(3:4)*q(3:4)' q2*q3-q1*q4 q2*q4+q1*q3;...
    q2*q3+q1*q4 -q([2 4])*q([2 4])' q3*q4-q1*q2;...
    q2*q4-q1*q3 q1*q2+q3*q4 -q(2:3)*q(2:3)'];
RRT = simplify(R*R');
Rq = simplify(jacobian(R(:),q));
Rdq = simplify(R'*reshape(Rq*qd',3,[]));
w = vex3(Rdq);
Wq = [-q2 q1 q4 -q3; -q3 -q4 q1 q2; -q4 q3 -q2 q1];
W = (q*q')/2*Wq'/(Wq*Wq');
W2 = -q*q'*[vex3(simplify(R'*reshape(Rq(:,1),3,3))) vex3(simplify(R'*reshape(Rq(:,2),3,3)))...
    vex3(simplify(R'*reshape(Rq(:,3),3,3))) vex3(simplify(R'*reshape(Rq(:,4),3,3)))];
Rq2 = reshape(Rq,3,3,[]);
% matlabFunction(Rq2,'File','quat_dR','Vars',{q'});
%%
clc; clear all;
q1 = rand(1,4); 
q2 = rand(1,4); 
q3 = quatmultiply(q1,q2);
q4 = quatx(q1',q2');
q5 = quatx21(q1',q2');
R1 = quat2rotm(q1);
R2 = quat2R(q1',q1*q1');
%% 
clc; clear all;
h = rand(4,1); v = rand(3,1);
nh2 = h'*h;
R = quat2R(h,nh2);
dR = quat_dR(h);
dRv1 = R*quat_dRv(h,nh2,v);
dRv2 = [dR(:,:,1)*v dR(:,:,2)*v...
    dR(:,:,3)*v dR(:,:,4)*v];
dRTv1 = quat_dRv(h,nh2,-(R'*v));
dRTv2 = [dR(:,:,1)'*v dR(:,:,2)'*v...
    dR(:,:,3)'*v dR(:,:,4)'*v];