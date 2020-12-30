function [mu,v]=gx2stat(lambda,m,delta,sigma,c)
	
	% GX2STAT Returns the mean and variance of a generalized chi-squared
	% variable (a weighted sum of non-central chi-squares).
	%
	% Abhranil Das <abhranil.das@utexas.edu>
	% Center for Perceptual Systems, University of Texas at Austin
	% If you use this code, please cite:
	% <a href="matlab:web('https://arxiv.org/abs/2012.14331')"
	% >A method to integrate and classify normal distributions</a>.
	%
	% Usage:
	% [mu,v]=gx2stat(lambda,m,delta,sigma,c)
	%
	% Example:
	% [mu,v]=gx2stat([1 -5 2],[1 2 3],[2 3 7],4,0)
	%
	% Required inputs:
	% lambda    row vector of coefficients of the non-central chi-squares
	% m         row vector of degrees of freedom of the non-central chi-squares
	% delta     row vector of non-centrality paramaters (sum of squares of
	%           means) of the non-central chi-squares
	% sigma     sd of normal term
	% c         constant term
	%
	% Outputs:
	% mu        mean
	% v         variance
    %
    % See also:
	% <a href="matlab:open(strcat(fileparts(which('gx2cdf')),filesep,'doc',filesep,'GettingStarted.mlx'))">Interactive demos</a>
	% gx2rnd
	% gx2_params_norm_quad
	
	parser = inputParser;
	addRequired(parser,'lambda',@(x) isreal(x) && isrow(x));
	addRequired(parser,'m',@(x) isreal(x) && isrow(x));
	addRequired(parser,'delta',@(x) isreal(x) && isrow(x));
	addRequired(parser,'c',@(x) isreal(x) && isscalar(x));
	parse(parser,lambda,m,delta,c);
	
	mu=dot(lambda,m+delta)+c;
	v=2*dot(lambda.^2,m+2*delta)+sigma^2;