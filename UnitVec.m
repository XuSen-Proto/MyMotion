function v = UnitVec(v)
v = v./(sqrt(sum(v.*v)));
end