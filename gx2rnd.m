function r=gx2rnd(w,k,lambda,s,m,varargin)

    % GX2RND Generates generalized chi-squared random numbers.
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
    % r=gx2rnd(w,k,lambda,s,m)
    % r=gx2rnd(w,k,lambda,s,m,sz)
    % r=gx2rnd(w,k,lambda,s,m,[sz1,sz2,...])
    % r=gx2rnd(w,k,lambda,s,m,sz,'method',method)
    %
    % Example:
    % r=gx2rnd([1 -5 2],[1 2 3],[2 3 7],5,1,5)
    %
    % Required inputs:
    % w         row vector of weights of the non-central chi-squares
    % k         row vector of degrees of freedom of the non-central chi-squares
    % lambda    row vector of non-centrality paramaters (sum of squares of
    %           means) of the non-central chi-squares
    % s         scale of normal term
    % m         offset
    %
    % Optional inputs:
    % sz        size(s) of the requested array
    % method    'sum' (default) to generate non-central chi-square and normal numbers and add them.
    %           'norm_quad' to generate standard normal vectors and compute their quadratic form.
    %
    % Output:
    % r         random number(s)
    %
    % See also:
    % <a href="matlab:open(strcat(fileparts(which('gx2cdf')),filesep,'doc',filesep,'GettingStarted.mlx'))">Getting Started guide</a>

    parser=inputParser;
    parser.KeepUnmatched=true;
    addRequired(parser,'w',@(x) isreal(x) && isrow(x));
    addRequired(parser,'k',@(x) isreal(x) && isrow(x));
    addRequired(parser,'lambda',@(x) isreal(x) && isrow(x));
    addRequired(parser,'s',@(x) isreal(x) && isscalar(x));
    addRequired(parser,'m',@(x) isreal(x) && isscalar(x));
    addOptional(parser,'sz',@(x) isreal(x));
    addParameter(parser,'method','sum',@(s) strcmpi(s,'sum') || strcmpi(s,'norm_quad'));
    parse(parser,w,k,lambda,s,m,varargin{:});
    sz=parser.Results.sz;
    if isscalar(sz), sz=[sz sz]; end
    method=parser.Results.method;

    if strcmpi(method,'sum')
        ncxs=arrayfun(@(w,k,lambda) w*ncx2rnd(k,lambda,sz),w,k,lambda,'un',0);

        r=zeros(size(ncxs{1}));
        for i=1:length(ncxs)
            r=r+ncxs{i};
        end
        r=r+normrnd(m,s,sz);
    elseif strcmpi(method,'norm_quad')
        % find the quadratic form of the standard multinormal
        quad=gx2_to_norm_quad_params(w,k,lambda,s,m);
        dim=numel(quad.q1);

        % generate standard normal vectors
        z=normrnd(0,1,[dim prod(sz,'all')]);

        % compute their quadratic form
        r=dot(z,quad.q2*z)+quad.q1'*z+quad.q0;
        r=reshape(r,sz);
    end

