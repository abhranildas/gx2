function [f,xgrid]=gx2pdf(x,w,k,lambda,m,s,varargin)

    % GX2PDF Returns the pdf of a generalized chi-squared (a weighted sum of
    % non-central chi-squares and a normal).
    %
    % Abhranil Das <abhranil.das@utexas.edu>
    % Center for Perceptual Systems, University of Texas at Austin
    % If you use this code, please cite:
    % <a href="matlab:web('https://arxiv.org/abs/2012.14331')"
    % >A method to integrate and classify normal distributions</a>.
    %
    % Usage:
    % f=gx2pdf(x,w,k,lambda,m,s)
    % f=gx2pdf(x,w,k,lambda,m,s,'dx',1e-3)
    % [f,xfull]=gx2pdf('full',w,k,lambda,m,s)
    % f=gx2pdf(x,w,k,lambda,m,s,'method','conv','AbsTol',0,'RelTol',1e-7)
    %
    % Example:
    % f=gx2pdf(25,[1 -5 2],[1 2 3],[2 3 7],5,0)
    % f=gx2pdf([17 25],[1 -5 2],[1 2 3],[2 3 7],5,0,'method','conv')
    % [f,xfull]=gx2pdf('full',[1 -5 2],[1 2 3],[2 3 7],5,0)
    %
    % Required inputs:
    % x         points at which to evaluate the pdf, or 'full', to return
    %           points xfull and f over a full range of the pdf (this uses
    %           the 'conv' method)
    % w         row vector of weights of the non-central chi-squares
    % k         row vector of degrees of freedom of the non-central chi-squares
    % lambda    row vector of non-centrality paramaters (sum of squares of
    %           means) of the non-central chi-squares
    % m         mean of normal term
    % s         sd of normal term
    %
    % Optional name-value inputs:
    % method    'diff' (default) for differentiating the generalized
    %           chi-square cdf, or 'conv' for convolving noncentral
    %           chi-square pdf's
    % dx        step-size of fineness (for convolving or differentiating)
    % AbsTol    absolute error tolerance for the output when using 'diff'
    % RelTol    relative error tolerance for the output when using 'diff'
    %           The absolute OR the relative tolerance is satisfied.
    %
    % Output:
    % f         computed pdf at points, or over full range
    % xfull     x-values over the full pdf range, returned when x is 'full'
    %
    % See also:
    % <a href="matlab:open(strcat(fileparts(which('gx2cdf')),filesep,'doc',filesep,'GettingStarted.mlx'))">Interactive demos</a>
    % gx2cdf
    % gx2cdf_davies
    % gx2cdf_imhof
    % gx2cdf_ruben

    parser = inputParser;
    addRequired(parser,'x',@(x) isreal(x) || strcmpi(x,'full'));
    addRequired(parser,'w',@(x) isreal(x) && isrow(x));
    addRequired(parser,'k',@(x) isreal(x) && isrow(x));
    addRequired(parser,'lambda',@(x) isreal(x) && isrow(x));
    addRequired(parser,'m',@(x) isreal(x) && isscalar(x));
    addRequired(parser,'s',@(x) isreal(x) && isscalar(x));
    addParameter(parser,'method','auto');
    addParameter(parser,'AbsTol',1e-10,@(x) isreal(x) && isscalar(x) && (x>=0));
    addParameter(parser,'RelTol',1e-6,@(x) isreal(x) && isscalar(x) && (x>=0));
    addParameter(parser,'xrange',[],@(x) isscalar(x) && (x>0));
    [~,v]=gx2stat(w,k,lambda,m,s);
    addParameter(parser,'n_grid',1e4+1,@(x) isscalar(x) && (x>0));
    addParameter(parser,'dx',sqrt(v)/1e4,@(x) isreal(x) && isscalar(x) && (x>=0)); % default derivative step-size is sd/100.
    parse(parser,x,w,k,lambda,m,s,varargin{:});

    method=parser.Results.method;

    if strcmpi(x,'full')
        method='conv';
    end

    if ~s && length(unique(w))==1 && ~strcmpi(x,'full')
        f=ncx2pdf((x-m)/unique(w),sum(k),sum(lambda))/abs(unique(w));
    elseif strcmpi(method,'ifft')
        [f,xgrid]=gx2pdf_ifft(x,w,k,lambda,m,s,varargin{:});
    elseif strcmpi(method,'conv')
        [f,xgrid]=gx2pdf_conv(x,w,k,lambda,m,s,varargin{:});
    elseif strcmpi(method,'diff')
        p_left=gx2cdf(x-dx,w,k,lambda,m,s,varargin{:});
        p_right=gx2cdf(x+dx,w,k,lambda,m,s,varargin{:});
        f=(p_right-p_left)/(2*dx);
        f=max(f,0);
    end
end