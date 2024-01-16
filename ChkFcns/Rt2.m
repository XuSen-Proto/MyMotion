function t = Rt2(rn,t,Sth,c0)
t = so3upoly([1 Sth c0],rn)*t;
end