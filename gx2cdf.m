function p=gx2cdf(x,lambda,m,delta,sigma,c,varargin)
	
	% GX2CDF Returns the cdf of a generalized chi-squared (a weighted sum of
	% non-central chi-squares and a normal), using Ruben's [1962] method,
	% Davies' [1973] method, or the native ncx2cdf, depending on the input.
	%
	% Abhranil Das <abhranil.das@utexas.edu>
	% Center for Perceptual Systems, University of Texas at Austin
	% If you use this code, please cite:
	% <a href="matlab:web('https://arxiv.org/abs/2012.14331')"
	% >A method to integrate and classify normal distributions</a>.
	%
	% Usage:
	% p=gx2cdf(x,lambda,m,delta,sigma,c)
	% p=gx2cdf(x,lambda,m,delta,sigma,c,'upper')
	% p=gx2cdf(x,lambda,m,delta,sigma,c,'AbsTol',0,'RelTol',1e-7)
	%
	% Example:
	% f=gx2pdf(25,[1 -5 2],[1 2 3],[2 3 7],5,0)
	%
	% Required inputs:
	% x         points at which to evaluate the CDF
	% lambda    row vector of coefficients of the non-central chi-squares
	% m         row vector of degrees of freedom of the non-central chi-squares
	% delta     row vector of non-centrality paramaters (sum of squares of
	%           means) of the non-central chi-squares
	% sigma     sd of normal term
	% c         constant term
	%
	% Optional positional input:
	% 'upper'   more accurate estimate of the complementary CDF when it's small
	%
	% Optional name-value inputs:
	% 'AbsTol'  absolute error tolerance for the output
	% 'RelTol'  relative error tolerance for the output
	%           The absolute OR the relative tolerance is satisfied.
	%
	% Output:
	% p         computed cdf
    %
	% See also:
	% <a href="matlab:open(strcat(fileparts(which('gx2cdf')),filesep,'doc',filesep,'GettingStarted.mlx'))">Interactive demos</a>
	% gx2cdf_davies
	% gx2cdf_imhof
    % gx2cdf_ruben
    % gx2pdf
	
	parser = inputParser;
	parser.KeepUnmatched = true;
	addRequired(parser,'x',@(x) isreal(x));
	addRequired(parser,'lambda',@(x) isreal(x) && isrow(x));
	addRequired(parser,'m',@(x) isreal(x) && isrow(x));
	addRequired(parser,'delta',@(x) isreal(x) && isrow(x));
	addRequired(parser,'sigma',@(x) isreal(x) && isscalar(x));
	addRequired(parser,'c',@(x) isreal(x) && isscalar(x));
	addOptional(parser,'side','lower',@(x) strcmpi(x,'lower') || strcmpi(x,'upper') );
	addParameter(parser,'AbsTol',1e-10,@(x) isreal(x) && isscalar(x) && (x>=0));
	addParameter(parser,'RelTol',1e-6,@(x) isreal(x) && isscalar(x) && (x>=0));
	
	parse(parser,x,lambda,m,delta,sigma,c,varargin{:});
	side=parser.Results.side;
	
	if ~sigma && length(unique(lambda))==1
		% native ncx2 fallback
		if (sign(unique(lambda))==1 && strcmpi(side,'lower')) || (sign(unique(lambda))==-1 && strcmpi(side,'upper'))
			p=ncx2cdf((x-c)/unique(lambda),sum(m),sum(delta));
		else
			p=ncx2cdf((x-c)/unique(lambda),sum(m),sum(delta),'upper');
		end
	elseif ~sigma && (all(lambda>0)||all(lambda<0))
		try
			p=gx2cdf_ruben(x,lambda,m,delta,c,varargin{:});
		catch
			p=gx2cdf_davies(x,lambda,m,delta,sigma,c,varargin{:});
		end
	else
		p=gx2cdf_davies(x,lambda,m,delta,sigma,c,varargin{:});
	end