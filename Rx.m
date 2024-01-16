function R = Rx(varargin)
if nargin==1
    th = varargin{1};
    sth = sin(th); cth = cos(th);
else
    cth = varargin{1};
    sth = varargin{2};
end
R = [1 0 0; 0 cth -sth; 0 sth cth];
end