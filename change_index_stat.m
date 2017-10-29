function ci = change_index_stat(a, b)

% check inputs
if isempty(a)
	error('%s: a is empty', mfilename)
elseif isempty(b)
	error('%s: a is empty', mfilename)
elseif length(a) ~= length(b)	
	error('%s: mismatch in length of a and b', mfilename);
elseif ~isnumeric(a) || ~isnumeric(b)
	error('%s: a and b must be numeric arrays', mfilename);
else
	% force to row vectors
	a = force_row(a);
	b = force_row(b);
end

% compute index
ci = (b - a) ./ (b + a);

