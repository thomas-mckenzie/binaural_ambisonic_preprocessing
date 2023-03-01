function y = rotate_vect(x, nr)
%
%  y = rotate_vect(x, nr)
%
% Right circular rotation of x by nr places, thus y(nr+1) = x(1)
% and y(nr) = x(length(x)).
%
% Bill Gardner
% Copyright 1995 MIT Media Lab. All rights reserved.
%
n = length(x);
if (nr <= 0)
	nr = n + nr;
end
y = zeros(size(x));
y(1 : nr) = x(n - nr + 1 : n);
y(nr + 1 : n) = x(1 : n - nr);
