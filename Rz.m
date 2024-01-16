function R = Rz(varargin)
if nargin==1
    th = varargin{1};
    sth = sin(th); cth = cos(th);
else
    cth = varargin{1};
    sth = varargin{2};
end
R = [cth -sth 0; sth cth 0; 0 0 1];
end