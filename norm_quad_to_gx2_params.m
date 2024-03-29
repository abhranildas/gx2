function [w,k,lambda,m,s]=norm_quad_to_gx2_params(mu,v,quad)
	
	% NORM_QUAD_TO_GX2_PARAMS A quadratic form of a normal vector is distributed
	% as a generalized chi-squared. This function takes the multinormal parameters
	% and the quadratic coeffs and returns the parameters of the generalized
	% chi-squared.
	%
	% Abhranil Das <abhranil.das@utexas.edu>
	% Center for Perceptual Systems, University of Texas at Austin
	% If you use this code, please cite:
	% <a href="matlab:web('https://arxiv.org/abs/2012.14331')"
	% >A method to integrate and classify normal distributions</a>.
	%
	% Example:
	% mu=[1;2]; % mean
	% v=[2 1; 1 3]; % covariance matrix
	% % Say q(x)=(x1+x2)^2-x1-1 = [x1;x2]'*[1 1; 1 1]*[x1;x2] + [-1;0]'*[x1;x2] - 1:
	% quad.q2=[1 1; 1 1];
	% quad.q1=[-1;0];
	% quad.q0=-1;
	%
	% [w,k,lambda,m,s]=gx2_params_norm_quad(mu,v,quad)
	%
	% Required inputs:
	% mu        column vector of normal mean
	% v         normal covariance matrix
	% quad      struct with quadratic form coefficients:
	%               q2      matrix of quadratic coefficients
	%               q1      column vector of linear coefficients
	%               q0      scalar constant
	%
	% Outputs:
	% w         row vector of weights of the non-central chi-squares
	% k         row vector of degrees of freedom of the non-central chi-squares
	% lambda    row vector of non-centrality paramaters (sum of squares of
	%           means) of the non-central chi-squares
	% m         mean of normal term
    % s         sd of normal term
    %
    % See also:
	% <a href="matlab:open(strcat(fileparts(which('gx2cdf')),filesep,'doc',filesep,'GettingStarted.mlx'))">Interactive demos</a>
	% gx2rnd
	% gx2stat
    % gx2cdf
    % gx2pdf
	
	% standardize the space
	q2s=0.5*(quad.q2+quad.q2'); % symmetrize q2
	q2=sqrtm(v)*q2s*sqrtm(v);
	q1=sqrtm(v)*(2*q2s*mu+quad.q1);
	q0=mu'*q2s*mu+quad.q1'*mu+quad.q0;
	
	[R,D]=eig(q2);
	d=diag(D)';
	b=(R'*q1)';
	
	[w,~,ic]=uniquetol(nonzeros(d)'); % unique non-zero eigenvalues
	k=accumarray(ic,1)'; % total dof of each eigenvalue
	lambda=arrayfun(@(x) sum((b(d==x)).^2),w)./(4*w.^2); % total non-centrality for each eigenvalue
    m=q0-dot(w,lambda);
    s=norm(b(~d));