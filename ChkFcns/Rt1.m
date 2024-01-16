function t = Rt1(rn,t,Sth,c0)
t = t+Ahat3B(rn,Sth*t+Ahat3B(rn,c0*t));
end