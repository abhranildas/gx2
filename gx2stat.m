function [mu,v]=gx2stat(w,k,lambda,m,s)
	
	% GX2STAT Returns the mean and variance of a generalized chi-squared
	% variable (a weighted sum of non-central chi-squares).
	%
	% Abhranil Das <abhranil.das@utexas.edu>
	% Center for Perceptual Systems, University of Texas at Austin
	% If you use this code, please cite:
	% <a href="matlab:web('https://arxiv.org/abs/2012.14331')"
	% >A method to integrate and classify normal distributions</a>.
	%
	% Example:
	% [mu,v]=gx2stat([1 -5 2],[1 2 3],[2 3 7],4,0)
	%
	% Required inputs:
	% w         row vector of weights of the non-central chi-squares
	% k         row vector of degrees of freedom of the non-central chi-squares
	% lambda    row vector of non-centrality paramaters (sum of squares of
	%           means) of the non-central chi-squares
	% m         mean of normal term
    % s         sd of normal term
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
	addRequired(parser,'w',@(x) isreal(x) && isrow(x));
	addRequired(parser,'k',@(x) isreal(x) && isrow(x));
	addRequired(parser,'lambda',@(x) isreal(x) && isrow(x));
	addRequired(parser,'m',@(x) isreal(x) && isscalar(x));
	addRequired(parser,'s',@(x) isreal(x) && isscalar(x));
	parse(parser,w,k,lambda,m,s);
	
	mu=dot(w,k+lambda)+m;
	v=2*dot(w.^2,k+2*lambda)+s^2;