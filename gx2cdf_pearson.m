function [p,flag]=gx2cdf_pearson(x,w,k,lambda,m,varargin)

    % GX2CDF_IMHOF Returns the cdf of a generalized chi-squared (a weighted sum of
    % non-central chi-squares), using Imhof's [1961] method.
    %
    % Abhranil Das <abhranil.das@utexas.edu>
    % Center for Perceptual Systems, University of Texas at Austin
    % If you use this code, please cite:
    % <a href="matlab:web('https://arxiv.org/abs/2012.14331')"
    % >A method to integrate and classify normal distributions</a>.
    %
    % Usage:
    % p=gx2cdf_imhof(x,w,k,lambda,m)
    % p=gx2cdf_imhof(x,w,k,lambda,m,'upper')
    % p=gx2cdf_imhof(x,w,k,lambda,m,'AbsTol',0,'RelTol',1e-7)
    % p=gx2cdf_imhof(x,w,k,lambda,m,'upper','approx','tail')
    %
    % Example:
    % p=gx2cdf_imhof(25,[1 -5 2],[1 2 3],[2 3 7],0)
    %
    % Required inputs:
    % x         points at which to evaluate the cdf
    % w         row vector of weights of the non-central chi-squares
    % k         row vector of degrees of freedom of the non-central chi-squares
    % lambda    row vector of non-centrality paramaters (sum of squares of
    %           means) of the non-central chi-squares
    % m         mean of normal term
    %
    % Optional positional input:
    % 'upper'   more accurate estimate of the complementary CDF when it's small
    %
    % Optional name-value inputs:
    % 'AbsTol'  absolute error tolerance for the output
    % 'RelTol'  relative error tolerance for the output
    %           The absolute OR the relative tolerance is satisfied.
    % 'approx'  set to 'tail' for Pearson's approximation of the tail
    %           probabilities. Works best for the upper (lower) tail when all
    %           w are positive (negative).
    %
    % Outputs:
    % p         computed cdf
    % flag      =true if output was too close to 0 or 1 to compute exactly with
    %           default settings. Try stricter tolerances or tail approx. for
    %           more accuracy.
    %
    % See also:
    % <a href="matlab:open(strcat(fileparts(which('gx2cdf')),filesep,'doc',filesep,'GettingStarted.mlx'))">Interactive demos</a>
    % gx2cdf_davies
    % gx2cdf_ruben
    % gx2cdf

    parser = inputParser;
    addRequired(parser,'x',@(x) isreal(x));
    addRequired(parser,'w',@(x) isreal(x) && isrow(x));
    addRequired(parser,'k',@(x) isreal(x) && isrow(x));
    addRequired(parser,'lambda',@(x) isreal(x) && isrow(x));
    addRequired(parser,'m',@(x) isreal(x) && isscalar(x));
    addOptional(parser,'side','lower',@(x) strcmpi(x,'lower') || strcmpi(x,'upper') );
    addParameter(parser,'output','cdf',@(x) strcmpi(x,'cdf') || strcmpi(x,'pdf') );

    parse(parser,x,w,k,lambda,m,varargin{:});
    side=parser.Results.side;

    j=(1:3)';
    c=sum((w.^j).*(j.*lambda+k),2);
    h=c(2)^3/c(3)^2;
    if c(3)>0
        y=(x-m-c(1))*sqrt(h/c(2))+h;
        if strcmpi(parser.Results.output,'cdf')
            if strcmpi(side,'lower')
                p=chi2cdf(y,h);
            elseif strcmpi(side,'upper')
                p=chi2cdf(y,h,'upper');
            end
        elseif strcmpi(parser.Results.output,'pdf')
            p=sqrt(h/c(2))*chi2pdf(y,h);
        end
    else
        c=sum(((-w).^j).*(j.*lambda+k),2);
        y=(-(x-m)-c(1))*sqrt(h/c(2))+h;
        if strcmpi(parser.Results.output,'cdf')
            if strcmpi(side,'lower')
                p=chi2cdf(y,h,'upper');
            elseif strcmpi(side,'upper')
                p=chi2cdf(y,h);
            end
        elseif strcmpi(parser.Results.output,'pdf')
            p=sqrt(h/c(2))*chi2pdf(y,h);
        end
    end


    flag = p<0 | p>1;
    p=max(p,0);
    p=min(p,1);

end