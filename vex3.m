% convert so3 matrix to R3 vector
function v = vex3(S)
v = 0.5*[S(3,2)-S(2,3); S(1,3)-S(3,1); S(2,1)-S(1,2)];
end