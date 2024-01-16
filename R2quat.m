function h = R2quat(R)
cthh = sqrt((1+R(1,1)+R(2,2)+R(3,3))/4);
if cthh~=0
    asthh = vex3(R-R')./(4*cthh);
else
    [~,im] = max([R(1,1) R(2,2) R(3,3)]);
    asthh = R(:,im);
    asthh(im) = asthh(im)+1;
    asthh = asthh./sqrt(2*asthh(im));
end
h = [cthh; asthh];
end