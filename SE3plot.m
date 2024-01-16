function t = SE3plot(T,FrameSize,ParentAx)
h = hggroup;
plot3(h,[0 FrameSize], [0 0], [0 0],'-r',...
    [0 0], [0 FrameSize], [0 0],'-g',...
    [0 0], [0 0], [0 FrameSize],'-b');
t = hgtransform('Parent',ParentAx);
set(h,'Parent',t);
set(t,'Matrix',T)
drawnow
end