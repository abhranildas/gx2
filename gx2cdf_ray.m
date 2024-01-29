function p=gx2cdf_ray(x,w,k,lambda,m,s,varargin)

    % GX2CDF_RAY Returns the cdf of a generalized chi-squared (a weighted
    % sum of non-central chi-squares and a normal), using the ray method [CITE].
    %
    % Abhranil Das <abhranil.das@utexas.edu>
    % Center for Perceptual Systems, University of Texas at Austin
    % If you use this code, please cite:
    % <a href="matlab:web('https://arxiv.org/abs/2012.14331')"
    % >A method to integrate and classify normal distributions</a>.
    %
    % Usage:
    % p=gx2cdf_davies(x,w,k,lambda,m,s)
    % p=gx2cdf_davies(x,w,k,lambda,m,s,'upper')
    % p=gx2cdf_davies(x,w,k,lambda,m,s,'AbsTol',0,'RelTol',1e-7)
    %
    % Example:
    % p=gx2cdf_davies(25,[1 -5 2],[1 2 3],[2 3 7],5,0)
    %
    % Required inputs:
    % x         points at which to evaluate the cdf
    % w         row vector of weights of the non-central chi-squares
    % k         row vector of degrees of freedom of the non-central chi-squares
    % lambda    row vector of non-centrality paramaters (sum of squares of
    %           means) of the non-central chi-squares
    % m         mean of normal term
    % s         sd of normal term
    %
    % Optional positional input:
    % 'upper'   more accurate estimate of the complementary CDF when it's small
    %
    % Optional name-value inputs:
    % 'AbsTol'  absolute error tolerance for the output
    % 'RelTol'  relative error tolerance for the output
    %           The absolute OR the relative tolerance is satisfied.
    %
    % Outputs:
    % p         computed cdf
    % flag      =true if output was too close to 0 or 1 to compute exactly with
    %           default settings. Try stricter tolerances.
    %
    % See also:
    % <a href="matlab:open(strcat(fileparts(which('gx2cdf')),filesep,'doc',filesep,'GettingStarted.mlx'))">Interactive demos</a>
    % gx2cdf_imhof
    % gx2cdf_ruben
    % gx2cdf

    parser=inputParser;
    parser.KeepUnmatched=true;
    addRequired(parser,'x',@(x) isreal(x));
    addRequired(parser,'w',@(x) isreal(x) && isrow(x));
    addRequired(parser,'k',@(x) isreal(x) && isrow(x));
    addRequired(parser,'lambda',@(x) isreal(x) && isrow(x));
    addRequired(parser,'m',@(x) isreal(x) && isscalar(x));
    addRequired(parser,'s',@(x) isreal(x) && isscalar(x));
    addOptional(parser,'side','lower',@(x) strcmpi(x,'lower') || strcmpi(x,'upper') );
    addParameter(parser,'AbsTol',1e-10,@(x) isreal(x) && isscalar(x) && (x>=0));
    addParameter(parser,'RelTol',1e-2,@(x) isreal(x) && isscalar(x) && (x>=0));
    addParameter(parser,'mc_samples',500);

    parse(parser,x,w,k,lambda,m,s,varargin{:});
    side=parser.Results.side;
    AbsTol=parser.Results.AbsTol;
    RelTol=parser.Results.RelTol;
    mc_samples=parser.Results.mc_samples;

    % find the quadratic form of the standard multinormal
    quad=gx2_to_norm_quad_params(w,k,lambda,m,s);
    dim=numel(quad.q1);

    % flip the domain, for lower side
    if strcmpi(side,'lower')
        quad=structfun(@uminus,quad,'un',0);
        x=-x;
    end

    % function to compute cdf for each point in x
    function p_each=gx2cdf_ray_each(x_each,quad,dim,AbsTol,RelTol,mc_samples)
        tic
        quad.q0=quad.q0-x_each;
        p_each=int_norm_ray(zeros(dim,1),eye(dim),quad,'AbsTol',AbsTol,'RelTol',RelTol,'mc_samples',mc_samples);
%         p_each=integrate_normal(zeros(dim,1),eye(dim),quad,'method','ray','AbsTol',AbsTol,'RelTol',RelTol,'mc_samples',mc_samples,'plotmode',0);
        toc
    end


    % integrate the standard multinormal over the domain using the ray method
    p=arrayfun(@(x_each) gx2cdf_ray_each(x_each,quad,dim,AbsTol,RelTol,mc_samples),x,'un',0);

    % if p is a cell but there are no symbols, convert to numeric array
    if iscell(p)
        num_idx=cellfun(@isnumeric, p);
        if all(num_idx)
            p=cell2mat(p);
        end
    end
end