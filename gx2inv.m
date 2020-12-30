function x=gx2inv(p,lambda,m,delta,sigma,c,varargin)
	
	% GX2INV Returns the inverse cdf of a generalized chi-squared, using
	% Ruben's [1962] method, Davies' [1973] method, or the native ncx2inv,
	% depending on the input.
	%
	% Abhranil Das <abhranil.das@utexas.edu>
	% Center for Perceptual Systems, University of Texas at Austin
	% If you use this code, please cite:
	% <a href="matlab:web('https://arxiv.org/abs/2012.14331')"
	% >A method to integrate and classify normal distributions</a>.
	%
	% Usage:
	% x=gx2inv(p,lambda,m,delta,sigma,c)
	% x=gx2inv(p,lambda,m,delta,sigma,c,'AbsTol',0,'RelTol',1e-7)
	%
	% Example:
	% x=gx2inv(0.9,[1 -5 2],[1 2 3],[2 3 7],5,0)
	%
	% Required inputs:
	% p         probabilities at which to evaluate the inverse cdf
	% lambda    row vector of coefficients of the non-central chi-squares
	% m         row vector of degrees of freedom of the non-central chi-squares
	% delta     row vector of non-centrality paramaters (sum of squares of
	%           means) of the non-central chi-squares
	% sigma     sd of normal term
	% c         constant term
	%
	% Optional name-value inputs:
	% 'AbsTol'  absolute error tolerance for the cdf function that is inverted
	% 'RelTol'  relative error tolerance for the cdf function that is inverted
	%           The absolute OR the relative tolerance is satisfied.
	%
	% Output:
	% x         computed point
    %
    % See also:
	% <a href="matlab:open(strcat(fileparts(which('gx2cdf')),filesep,'doc',filesep,'GettingStarted.mlx'))">Interactive demos</a>
	% gx2cdf
    % gx2pdf
	% gx2rnd
    % gx2stat
	
	parser=inputParser;
	parser.KeepUnmatched=true;
	addRequired(parser,'p',@(x) isreal(x) && all(x>=0 & x<=1));
	addRequired(parser,'lambda',@(x) isreal(x) && isrow(x));
	addRequired(parser,'m',@(x) isreal(x) && isrow(x));
	addRequired(parser,'delta',@(x) isreal(x) && isrow(x));
	addRequired(parser,'sigma',@(x) isreal(x) && isscalar(x));
	addRequired(parser,'c',@(x) isreal(x) && isscalar(x));
	
	parse(parser,p,lambda,m,delta,sigma,c,varargin{:});
	
	if ~sigma && length(unique(lambda))==1
		% native ncx2 fallback
		if sign(unique(lambda))==1
			x=ncx2inv(p,sum(m),sum(delta))*unique(lambda)+c;
		elseif sign(unique(lambda))==-1
			x=ncx2inv(1-p,sum(m),sum(delta))*unique(lambda)+c;
		end
	else
		mu=gx2stat(lambda,m,delta,sigma,c);
		x=arrayfun(@(p) fzero(@(x) gx2cdf(x,lambda,m,delta,sigma,c,varargin{:})-p,mu),p);
	end