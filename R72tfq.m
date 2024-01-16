%convert a R7 vector into a quat x R3 element
function tf = R72tfq(tf)
tf = [UnitVec(tf(1:4,:)); tf(5:7,:)];
end