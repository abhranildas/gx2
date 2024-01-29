function f=gx2pdf_ray_fast(x,w,k,lambda,m,s,varargin)

    % GX2CDF_DAVIES Returns the cdf of a generalized chi-squared (a weighted
    % sum of non-central chi-squares and a normal), using Davies' [1973]
    % method, which is an extension of Imhof's [1961] method to include the
    % s term, so Imhof's method is not separately available in the toolbox
    % any more.
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
    % flag      = true for output(s) too close to 0 or 1 to compute exactly with
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
    addParameter(parser,'mc_samples',1e4);

    parse(parser,x,w,k,lambda,m,s,varargin{:});
    mc_samples=parser.Results.mc_samples;

    y=x(:); % flatten input array x into a column vector of levels y

    % find standard quadratic form corresponding to the gx2:
    quad=gx2_to_norm_quad_params(w,k,lambda,m,s);
    dim=numel(quad.q1);

    % uniform random rays (points on n-sphere)
    n_z=mvnrnd(zeros(dim,1),eye(dim),mc_samples)';
    n_z=n_z./vecnorm(n_z,2);

    % find the quadratic coefficients across all rays
    q2=dot(n_z,quad.q2*n_z);
    q1=quad.q1'*n_z;
    q0=quad.q0;

    % valid levels y for which there are roots
    valid_idx=(sign(q2).*y > sign(q2).*(q0-q1.^2./(4*q2)));

    % discriminant of the quadratic across rays and levels
    delta2=q1.^2-4*q2.*(q0-y); % delta^2
    delta=nan(size(delta2));
    delta(valid_idx)=sqrt(delta2(valid_idx));

    % roots across rays
    z=(-q1+cat(3,-1,1).*delta)./(2*q2);

    % 1. phi_ray at all roots across all rays
    % 2. sum across two roots on each ray
    % 3. divide by quadratic slope
    % 4. add across rays, gives us pdf at each level
    f=sum(sum(phi_ray(z,dim),3)./delta,2,'omitnan')/mc_samples;

    % reshape flattened array to shape of input x
    f=reshape(f,size(x));
end