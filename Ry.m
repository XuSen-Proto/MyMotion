function R = Ry(varargin)
if nargin==1
    th = varargin{1};
    sth = sin(th); cth = cos(th);
else
    cth = varargin{1};
    sth = varargin{2};
end
R = [cth 0 sth; 0 1 0; -sth 0 cth];
end