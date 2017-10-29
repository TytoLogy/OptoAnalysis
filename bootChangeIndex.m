function varargout = bootChangeIndex(A, B, nIter)

if isempty(A) && isempty(B)
	warning('A, B are empty, using diagnostir dummy data')
	
	nA = query_uservalue('# of test values', [5 1000 10]);
	nB = nA;
	muA = query_uservalue('Mean value, test data A', [-1000 1000 0]);
	sigmaA = query_uservalue('Std. Dev. value, test data A', [-1000 1000 1]);
	muB = query_uservalue('Mean value, test data B', [-1000 1000 10]);
	sigmaB = query_uservalue('Std. Dev. value, test data B', [-1000 1000 1]);
	%  generate test data
	A = normrnd(muA, sigmaA, [nA, 1]);
	B = normrnd(muB, sigmaB, [nB, 1]);
elseif isempty(A) || isempty(B)
	error('%s: A or B is empty!', mfilename);
end
% check nIter value
if isempty(nIter) || (nIter < 1)
	error('%s: invalid nIter! (empty or < 1)');
elseif nIter > 1e9
	error('%s: invalid nIter (way too big)!');
end

% some basic stats on inputs
statsIn.meanA = mean(A);
statsIn.stdA = std(A);
statsIn.medianA = median(A);
statsIn.meanB = mean(B);
statsIn.stdB = std(B);
statsIn.medianB = median(B);

% get samples from each sample
sampleA = datasample(A, nIter);
sampleB = datasample(B, nIter);

% compute change_index
sampleC = change_index_stat(sampleA, sampleB);

% compute output stats
statsOut.meanA = mean(sampleA);
statsOut.stdA = std(sampleA);
statsOut.medianA = median(sampleA);
statsOut.meanB = mean(sampleB);
statsOut.stdB = std(sampleB);
statsOut.medianB = median(sampleB);
statsOut.meanC = mean(sampleC);
statsOut.stdC = std(sampleC);
statsOut.medianC = median(sampleC);

% build output
out.inputstats = statsIn;
out.outputstats = statsOut;
out.sampleA = sampleA;
out.sampleB = sampleB;
out.sampleC = sampleC;

varargout{1} = out;
