function quad=gx2_to_norm_quad_params(w,k,lambda,s,m)

    % GX2_TO_NORM_QUAD_PARAMS A generalized chi-square variable is
    % a quadratic form of a standard normal vector. This function takes
    % the parameters of the generalized chi-square distribution, and finds
    % the coefficients of such a quadratic form.
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
    % Usage:
    % quad=gx2_to_norm_quad_params(w,k,lambda,s,m)
    %
    % Example:
    % quad=gx2_to_norm_quad_params([1 -5 2],[1 2 3],[2 3 7],10,5);
    %
    % Required inputs:
    % w         row vector of weights of the non-central chi-squares
    % k         row vector of degrees of freedom of the non-central chi-squares
    % lambda    row vector of non-centrality paramaters (sum of squares of
    %           means) of the non-central chi-squares
    % s         scale of normal term
    % m         offset
    %
    %
    % Outputs:
    % quad      struct with quadratic form coefficients:
    %               q2      matrix of quadratic coefficients
    %               q1      column vector of linear coefficients
    %               q0      scalar constant
    %           The dimension of the standard normal is the length of q1.
    %
    % See also:
    % <a href="matlab:open(strcat(fileparts(which('gx2cdf')),filesep,'doc',filesep,'GettingStarted.mlx'))">Getting Started guide</a>

    q2_parts=arrayfun(@(w,k) w*ones(1,k),w,k,'un',0); % each w_i, k_i times
    q2=horzcat(q2_parts{:}); % append all of them

    % each w_i*sqrt(lambda_i), append 0 k_i-1 times
    q1_parts=arrayfun(@(w,lambda,k) [w*sqrt(lambda) zeros(1,k-1)],w,lambda,k,'un',0);
    q1=-2*horzcat(q1_parts{:})'; % append all of them, multiply by -2

    if s % if there is a normal term,
        q2=[q2 0];
        q1=[q1;s];
    end

    quad.q2=diag(q2);
    quad.q1=q1;
    quad.q0=dot(w,lambda)+m;

