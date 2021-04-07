function r=gx2rnd(w,k,lambda,m,s,varargin)
	
	% GX2RND Generates generalized chi-squared  random numbers.
	%
	% Abhranil Das <abhranil.das@utexas.edu>
	% Center for Perceptual Systems, University of Texas at Austin
	% If you use this code, please cite:
	% <a href="matlab:web('https://arxiv.org/abs/2012.14331')"
	% >A method to integrate and classify normal distributions</a>.
	%
	% Usage:
	% r=gx2rnd(w,k,lambda,m,s)
	% r=gx2rnd(w,k,lambda,m,s,sz)
	% r=gx2rnd(w,k,lambda,m,s,sz1,sz2,...)
	% r=gx2rnd(w,k,lambda,m,s,[sz1,sz2,...])
	%
	% Example:
	% r=gx2rnd([1 -5 2],[1 2 3],[2 3 7],5,1,5)
	%
	% Required inputs:
	% w         row vector of weights of the non-central chi-squares
	% k         row vector of degrees of freedom of the non-central chi-squares
	% lambda    row vector of non-centrality paramaters (sum of squares of
	%           means) of the non-central chi-squares
	% m         mean of normal term
    % s         sd of normal term
	%
	% Optional positional input:
	% sz        size(s) of the requested array
	%
	% Output:
	% r         random number(s)
    %
    % See also:
	% <a href="matlab:open(strcat(fileparts(which('gx2cdf')),filesep,'doc',filesep,'GettingStarted.mlx'))">Interactive demos</a>
	% gx2stat
	% gx2_params_norm_quad
	
	parser=inputParser;
	parser.KeepUnmatched=true;
	addRequired(parser,'w',@(x) isreal(x) && isrow(x));
	addRequired(parser,'k',@(x) isreal(x) && isrow(x));
	addRequired(parser,'lambda',@(x) isreal(x) && isrow(x));
	addRequired(parser,'m',@(x) isreal(x) && isscalar(x));
	addRequired(parser,'s',@(x) isreal(x) && isscalar(x));
	parse(parser,w,k,lambda,m,s);
	
	ncxs=arrayfun(@(l,k,d) l*ncx2rnd(k,d,varargin{:}),w,k,lambda,'un',0);
	
	r=zeros(size(ncxs{1}));
	for i=1:length(ncxs)
		r=r+ncxs{i};
	end
	r=r+normrnd(m,s,varargin{:});
	
