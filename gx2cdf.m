function p=gx2cdf(x,w,k,lambda,m,s,varargin)

    % GX2CDF Returns the cdf of a generalized chi-squared (a weighted sum of
    % non-central chi-squares and a normal), using Ruben's [1962] method,
    % Davies' [1973] method, or the native ncx2cdf, depending on the input.
    %
    % Abhranil Das <abhranil.das@utexas.edu>
    % Center for Perceptual Systems, University of Texas at Austin
    % If you use this code, please cite:
    % <a href="matlab:web('https://arxiv.org/abs/2012.14331')"
    % >A method to integrate and classify normal distributions</a>.
    %
    % Usage:
    % p=gx2cdf(x,w,k,lambda,m,s)
    % p=gx2cdf(x,w,k,lambda,m,s,'upper')
    % p=gx2cdf(x,w,k,lambda,m,s,'AbsTol',0,'RelTol',1e-7)
    %
    % Example:
    % f=gx2cdf(25,[1 -5 2],[1 2 3],[2 3 7],5,0)
    %
    % Required inputs:
    % x         points at which to evaluate the CDF
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
    % AbsTol    absolute error tolerance for the output
    % RelTol    relative error tolerance for the output
    %           The absolute OR the relative tolerance is satisfied.
    % vpa       false (default) to do the integral in Imhof's method numerically,
    %           true to do it symbolically with variable precision.
    %
    % Output:
    % p         computed cdf
    %
    % See also:
    % <a href="matlab:open(strcat(fileparts(which('gx2cdf')),filesep,'doc',filesep,'GettingStarted.mlx'))">Interactive demos</a>
    % gx2cdf_davies
    % gx2cdf_imhof
    % gx2cdf_ruben
    % gx2pdf

    parser = inputParser;
    parser.KeepUnmatched = true;
    addRequired(parser,'x',@(x) isreal(x));
    addRequired(parser,'w',@(x) isreal(x) && isrow(x));
    addRequired(parser,'k',@(x) isreal(x) && isrow(x));
    addRequired(parser,'lambda',@(x) isreal(x) && isrow(x));
    addRequired(parser,'m',@(x) isreal(x) && isscalar(x));
    addRequired(parser,'s',@(x) isreal(x) && isscalar(x));
    addOptional(parser,'side','lower',@(x) strcmpi(x,'lower') || strcmpi(x,'upper') );
    addParameter(parser,'method','auto');
    addParameter(parser,'AbsTol',1e-10,@(x) isreal(x) && isscalar(x) && (x>=0));
    addParameter(parser,'RelTol',1e-6,@(x) isreal(x) && isscalar(x) && (x>=0));

    parse(parser,x,w,k,lambda,m,s,varargin{:});
    method=parser.Results.method;
    side=parser.Results.side;

    if strcmpi(method,'auto')
        if ~s && length(unique(w))==1
            % ncx2 fallback
            if (sign(unique(w))==1 && strcmpi(side,'lower')) || (sign(unique(w))==-1 && strcmpi(side,'upper'))
                p=ncx2cdf((x-m)/unique(w),sum(k),sum(lambda));
            else
                p=ncx2cdf((x-m)/unique(w),sum(k),sum(lambda),'upper');
            end
        else
            p=gx2cdf_ray(x,w,k,lambda,m,s,varargin{:});
        end
    elseif strcmpi(method,'ray')
        p=gx2cdf_ray(x,w,k,lambda,m,s,varargin{:});
    elseif strcmpi(method,'imhof')
        p=gx2cdf_imhof(x,w,k,lambda,m,s,varargin{:});
    elseif strcmpi(method,'ruben')
        if s || ~(all(w>0)||all(w<0))
            error("Ruben's method can only be used when all w are the same sign and s=0.")
        else
            p=gx2cdf_ruben(x,w,k,lambda,m,varargin{:});
        end
    end