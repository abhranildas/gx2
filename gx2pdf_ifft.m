function [f,xgrid]=gx2pdf_ifft(x,w,k,lambda,m,s,varargin)

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
    % range of x over which to compute ifft.
    % default span is upto 4 sd from mean
    [mu,v]=gx2stat(w,k,lambda,m,s);
    addParameter(parser,'xrange',max(abs(mu+[-1 1]*4*sqrt(v))),@(x) isscalar(x) && (x>0));
    % number of grid points for ifft
    addParameter(parser,'n_grid',1e4+1,@(x) isscalar(x) && (x>0));
    parse(parser,x,w,k,lambda,m,s,varargin{:});
    
    n_grid=round(parser.Results.n_grid);
    % make n_grid odd
    if mod(n_grid,2)
        n_grid=n_grid+1;
    end

    n=(n_grid-1)/2;
    idx=-n:n;

    xrange=parser.Results.xrange;
    dx=xrange/n;

    dt=2*pi/(n_grid*dx);

%     dt=.01;
    t=idx*dt; % grid of points to evaluate characteristic function

    phi=gx2char(-t,w,k,lambda,m,s);

    f=fftshift(ifft(ifftshift(phi)))*n_grid*dt/(2*pi);

    % grid of x over which this pdf is computed
%     dx=2*pi/(n_grid*dt);
    xgrid=idx*dx;

    % normalize
    f=f/(sum(f)*dx);

    if ~strcmpi(x,'full')
        F=griddedInterpolant(xgrid,f);
        f=F(x);
    end

    f=max(f,0);
end