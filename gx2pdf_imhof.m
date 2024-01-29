function [f,flag]=gx2pdf_imhof(x,w,k,lambda,m,s,varargin)

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
    addOptional(parser,'side','lower',@(x) strcmpi(x,'lower') || strcmpi(x,'upper') );
    addParameter(parser,'AbsTol',1e-10,@(x) isreal(x) && isscalar(x) && (x>=0));
    addParameter(parser,'RelTol',1e-6,@(x) isreal(x) && isscalar(x) && (x>=0));

    parse(parser,x,w,k,lambda,m,s,varargin{:});
    AbsTol=parser.Results.AbsTol;
    RelTol=parser.Results.RelTol;

    % compute the integral
    if any(strcmpi(parser.UsingDefaults,'AbsTol')) && any(strcmpi(parser.UsingDefaults,'RelTol'))
        imhof_integral=arrayfun(@(x) integral(@(u) gx2pdf_imhof_integrand(u,x,w',k',lambda',m,s),0,inf),x);
        f=imhof_integral/(2*pi);
    else
        syms u
        imhof_integral=arrayfun(@(x) vpaintegral(@(u) gx2pdf_imhof_integrand(u,x,w',k',lambda',m,s),...
            u,0,inf,'AbsTol',AbsTol,'RelTol',RelTol,'MaxFunctionCalls',inf),x);
        f=double(imhof_integral/(2*pi));
    end

    flag = f<0;
    if any(flag)
        warning('gx2pdf output(s) too close to 0 to compute exactly, so clipping to 0. Check the flag output, and try stricter tolerances.')
        f=max(f,0);
    end


end