function out = RandInt(varargin);
sigStr = '';
m = [];
n = [];
range = [];
state = [];
for i=1:nargin
if(i>1)
sigStr(size(sigStr,2)+1) = '/';
end;
if(isnumeric(varargin{i}))
sigStr(size(sigStr,2)+1) = 'n';
else
error('Only numeric arguments are accepted.');
end;
end;
switch sigStr
case ''
case 'n'
m		= varargin{1};
case 'n/n'
m		= varargin{1};
n		= varargin{2};
case 'n/n/n'
m		= varargin{1};
n  	= varargin{2};
range = varargin{3};
case 'n/n/n/n'
m		= varargin{1};
n		= varargin{2};
range = varargin{3};
state = varargin{4};
otherwise
error('Syntax error.');
end;
if isempty(m)
m = 1;
end
if isempty(n)
n = m;
end
if isempty(range)
range = [0, 1];
end
len_range = size(range,1) * size(range,2);
if (~isfinite(m)) | (~isfinite(n))
error('Matrix dimensions must be finite.');
elseif (floor(m) ~= m) | (floor(n) ~= n) | (~isreal(m)) | (~isreal(n))
error('Matrix dimensions must be real integers.');
elseif (m < 0) | (n < 0)
error('Matrix dimensions must be positive.');
elseif (length(m) > 1) | (length(n) > 1)
error('Matrix dimensions must be scalars.');
elseif len_range > 2
error('The RANGE parameter should contain no more than two elements.');
elseif max(max(floor(range) ~= range)) | (~isreal(range)) | (~isfinite(range))
error('The RANGE parameter must only contain real finite integers.');
end
if len_range < 2
if range < 0
range = [range+1, 0];
elseif range > 0
range = [0, range-1];
else
range = [0, 0];    
end
end
range = sort(range);
distance = range(2) - range(1);
if ~isempty(state)
rand('state', state);
end
r = floor(rand(m, n) * (distance+1));
out = ones(m,n)*range(1);
out = out + r;
