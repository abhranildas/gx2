function [p,errbnd]=gx2cdf_ruben(x,w,k,lambda,m,varargin)
	
	% GX2CDF_RUBEN Returns the cdf of a generalized chi-squared (a weighted sum of
	% non-central chi-squares with all weights the same sign), using Ruben's
	% [1962] method.
	%
	% Abhranil Das <abhranil.das@utexas.edu>
	% Center for Perceptual Systems, University of Texas at Austin
	% If you use this code, please cite:
	% <a href="matlab:web('https://arxiv.org/abs/2012.14331')"
	% >A method to integrate and classify normal distributions</a>.
	%
	% Usage:
	% p=gx2cdf_ruben(x,w,k,lambda,m)
	% p=gx2cdf_ruben(x,w,k,lambda,m,N)
	% [p,err]=gx2cdf_ruben(x,w,k,lambda,m)
	%
	% Example:
	% [p,err]=gx2cdf_ruben(25,[1 5 2],[1 2 3],[2 3 7],0,100)
	%
	% Required inputs:
	% x         points at which to evaluate the cdf
	% w         row vector of weights of the non-central chi-squares
	% k         row vector of degrees of freedom of the non-central chi-squares
	% lambda    row vector of non-centrality paramaters (sum of squares of
	%           means) of the non-central chi-squares
	% m         mean of normal term
	%
	% Optional positional input:
	% N         no. of terms in the approximation. Default = 1000.
	%
	% Outputs:
	% p         computed cdf
	% errbnd    upper error bound of the computed cdf
	%
    % See also:
	% <a href="matlab:open(strcat(fileparts(which('gx2cdf')),filesep,'doc',filesep,'GettingStarted.mlx'))">Interactive demos</a>
	% gx2cdf_davies
	% gx2cdf_imhof
    % gx2cdf

	parser=inputParser;
	parser.KeepUnmatched=true;
	addRequired(parser,'x',@(x) isreal(x));
	addRequired(parser,'w',@(x) isreal(x) && isrow(x)  && (all(x>0)||all(x<0)) );
	addRequired(parser,'k',@(x) isreal(x) && isrow(x));
	addRequired(parser,'lambda',@(x) isreal(x) && isrow(x));
	addRequired(parser,'m',@(x) isreal(x) && isscalar(x));
	addOptional(parser,'side','lower',@(x) strcmpi(x,'lower') || strcmpi(x,'upper') );
	addParameter(parser,'N',1e2,@(x) ismember(x,1:x));
	
	parse(parser,x,w,k,lambda,m,varargin{:});
	side=parser.Results.side;
	N=parser.Results.N;
	lambda_pos=true;
	
	if all(w<0)
		w=-w; x=-x; m=-m; lambda_pos=false;
	end
	beta=0.90625*min(w);
	M=sum(k);
	n=(1:N-1)';
	
	% compute the g's
	g=sum(k.*(1-beta./w).^n,2)+ beta*n.*((1-beta./w).^(n-1))*(lambda./w)';
	
	% compute the expansion coefficients
	a=nan(N,1);
	a(1)=sqrt(exp(-sum(lambda))*beta^M*prod(w.^(-k)));
	if a(1)<realmin
		error('Underflow error: some series coefficients are smaller than machine precision.')
	end
	for j=1:N-1
		a(j+1)=dot(flip(g(1:j)),a(1:j))/(2*j);
	end
	
	% compute the central chi-squared integrals
	[xg,mg]=meshgrid((x-m)/beta,M:2:M+2*(N-1));
	F=arrayfun(@(x,k) chi2cdf(x,k),xg,mg);
	
	% compute the integral
	p=a'*F;
	
	% flip if necessary
	if (lambda_pos && strcmpi(side,'upper')) || (~lambda_pos && strcmpi(side,'lower'))
		p=1-p;
	end
	
	% compute the truncation error
	errbnd=(1-sum(a))*chi2cdf((x-m)/beta,M+2*N);