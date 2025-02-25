function [w,k,lambda,s,m]=norm_quad_to_gx2_params(mu,v,quad)
	
	% NORM_QUAD_TO_GX2_PARAMS A quadratic form of a normal vector is distributed
	% as a generalized chi-squared. This function takes the multinormal parameters
	% and the quadratic coefficients and returns the parameters of the generalized
	% chi-squared.
	%
    % Abhranil Das
    % Center for Perceptual Systems, University of Texas at Austin
    % Comments, questions, bugs to abhranil.das@utexas.edu
    % If you use this code, please cite:
    % 1. <a href="matlab:web('https://arxiv.org/abs/2012.14331')"
    % >A method to integrate and classify normal distributions</a>
    % 2. <a href="matlab:web('https://arxiv.org/abs/2404.05062')"
    % >New methods for computing the generalized chi-square distribution</a>
	%
	% Example:
	% mu=[1;2]; % mean
	% v=[2 1; 1 3]; % covariance matrix
	% % Say q(x)=(x1+x2)^2-x1-1 = [x1;x2]'*[1 1; 1 1]*[x1;x2] + [-1;0]'*[x1;x2] - 1:
	% quad.q2=[1 1; 1 1];
	% quad.q1=[-1;0];
	% quad.q0=-1;
	%
	% [w,k,lambda,s,m]=gx2_params_norm_quad(mu,v,quad)
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
    % s         scale of normal term
    % m         offset
    %
    % See also:
    % <a href="matlab:open(strcat(fileparts(which('gx2cdf')),filesep,'doc',filesep,'GettingStarted.mlx'))">Getting Started guide</a>
	
	% standardize the space
	q2_sym=0.5*(quad.q2+quad.q2'); % symmetrize q2

    % compute sqrtm of v while avoiding small negative eigenvalues
    [R, D] = eig(v);
    d = diag(D);
    d(d < 0) = 0;              % Threshold any small negatives to zero
    sqrt_d = sqrt(d);          % Compute the square roots
    sqrt_v = R * diag(sqrt_d) * R';  % Reassemble the square root matrix

	q2=sqrt_v*q2_sym*sqrt_v;
    q2=(q2+q2')/2; % symmetrize q2 again
	q1=sqrt_v*(2*q2_sym*mu+quad.q1);
	q0=mu'*q2_sym*mu+quad.q1'*mu+quad.q0;
	
	[R,D]=eig(q2);
	d=diag(D)';
	b=(R'*q1)';
	
	[w,~,ic]=uniquetol(nonzeros(d)'); % unique non-zero eigenvalues
	k=accumarray(ic,1)'; % total dof of each eigenvalue
	lambda=arrayfun(@(x) sum((b(d==x)).^2),w)./(4*w.^2); % total non-centrality for each eigenvalue
    m=q0-dot(w,lambda);
    s=norm(b(~d));