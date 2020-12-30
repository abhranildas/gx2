function r=gx2rnd(lambda,m,delta,sigma,c,varargin)
	
	% GX2RND Generates generalized chi-squared  random numbers.
	%
	% Abhranil Das <abhranil.das@utexas.edu>
	% Center for Perceptual Systems, University of Texas at Austin
	% If you use this code, please cite:
	% <a href="matlab:web('https://arxiv.org/abs/2012.14331')"
	% >A method to integrate and classify normal distributions</a>.
	%
	% Usage:
	% r=gx2rnd(lambda,m,delta,sigma,c)
	% r=gx2rnd(lambda,m,delta,sigma,c,sz)
	% r=gx2rnd(lambda,m,delta,sigma,c,sz1,sz2,...)
	% r=gx2rnd(lambda,m,delta,sigma,c,[sz1,sz2,...])
	%
	% Example:
	% r=gx2rnd([1 -5 2],[1 2 3],[2 3 7],5,1,5)
	%
	% Required inputs:
	% lambda    row vector of coefficients of the non-central chi-squares
	% m         row vector of degrees of freedom of the non-central chi-squares
	% delta     row vector of non-centrality paramaters (sum of squares of
	%           means) of the non-central chi-squares
	% sigma     sd of normal term
	% c         constant term
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
	addRequired(parser,'lambda',@(x) isreal(x) && isrow(x));
	addRequired(parser,'m',@(x) isreal(x) && isrow(x));
	addRequired(parser,'delta',@(x) isreal(x) && isrow(x));
	addRequired(parser,'sigma',@(x) isreal(x) && isscalar(x));
	addRequired(parser,'c',@(x) isreal(x) && isscalar(x));
	parse(parser,lambda,m,delta,sigma,c);
	
	ncxs=arrayfun(@(l,m,d) l*ncx2rnd(m,d,varargin{:}),lambda,m,delta,'un',0);
	
	r=zeros(size(ncxs{1}));
	for i=1:length(ncxs)
		r=r+ncxs{i};
	end
	r=r+normrnd(c,sigma,varargin{:});
	
