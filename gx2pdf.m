function f=gx2pdf(x,lambda,m,delta,sigma,c,varargin)
	
	% GX2PDF Returns the pdf of a generalized chi-squared (a weighted sum of
	% non-central chi-squares and a normal).
	%
	% Abhranil Das <abhranil.das@utexas.edu>
	% Center for Perceptual Systems, University of Texas at Austin
	% If you use this code, please cite:
	% <a href="matlab:web('https://arxiv.org/abs/2012.14331')"
	% >A method to integrate and classify normal distributions</a>.
	%
	% Usage:
	% f=gx2pdf(x,lambda,m,delta,sigma,c)
	% f=gx2pdf(x,lambda,m,delta,sigma,c,'dx',1e-3)
	% f=gx2pdf(x,lambda,m,delta,sigma,c,'AbsTol',0,'RelTol',1e-7)
	%
	% Example:
	% f=gx2pdf(25,[1 -5 2],[1 2 3],[2 3 7],5,0)
	%
	% Required inputs:
	% x         points at which to evaluate the pdf
	% lambda    row vector of coefficients of the non-central chi-squares
	% m         row vector of degrees of freedom of the non-central chi-squares
	% delta     row vector of non-centrality paramaters (sum of squares of
	%           means) of the non-central chi-squares
	% sigma     sd of normal term
	% c         constant term
	%
	% Optional name-value inputs:
	% dx        step-size for numerically differentiating cdf
	% 'AbsTol'  absolute error tolerance for the output
	% 'RelTol'  relative error tolerance for the output
	%           The absolute OR the relative tolerance is satisfied.
	%
	% Output:
	% f         computed pdf
    % See also:
	% <a href="matlab:open(strcat(fileparts(which('gx2cdf')),filesep,'doc',filesep,'GettingStarted.mlx'))">Interactive demos</a>
	% gx2cdf
    % gx2cdf_davies
	% gx2cdf_imhof
    % gx2cdf_ruben
	
	parser = inputParser;
	addRequired(parser,'x',@(x) isreal(x));
	addRequired(parser,'lambda',@(x) isreal(x) && isrow(x));
	addRequired(parser,'m',@(x) isreal(x) && isrow(x));
	addRequired(parser,'delta',@(x) isreal(x) && isrow(x));
	addRequired(parser,'sigma',@(x) isreal(x) && isscalar(x));
	addRequired(parser,'c',@(x) isreal(x) && isscalar(x));
	addParameter(parser,'AbsTol',1e-10,@(x) isreal(x) && isscalar(x) && (x>=0));
	addParameter(parser,'RelTol',1e-6,@(x) isreal(x) && isscalar(x) && (x>=0));
	[~,v]=gx2stat(lambda,m,delta,sigma,c);
	addParameter(parser,'dx',sqrt(v)/100,@(x) isreal(x) && isscalar(x) && (x>=0)); % default derivative step-size is sd/100.
	parse(parser,x,lambda,m,delta,sigma,c,varargin{:});
	
	if ~sigma && length(unique(lambda))==1
		f=ncx2pdf((x-c)/unique(lambda),sum(m),sum(delta))/abs(unique(lambda));
	else
		dx=parser.Results.dx;
		p_left=gx2cdf(x-dx,lambda,m,delta,sigma,c,varargin{:});
		p_right=gx2cdf(x+dx,lambda,m,delta,sigma,c,varargin{:});
		f=max((p_right-p_left)/(2*dx),0);
	end
	
end