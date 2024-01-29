function [f,xgrid]=gx2pdf_conv(x,w,k,lambda,m,s,varargin)

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

    xrange=parser.Results.xrange;
        
    w_idx=find(w); % non-zero w indices

    if isempty(xrange) % no xrange provided

        % find the xrange over which to convolve
        % xrange is a length that covers the
        % widest of the constituent pdfs
        xrange=nan(1,nnz(w));
        for i=1:nnz(w)
            xrange(i)=gx2inv(1-eps,abs(w(w_idx(i))),k(w_idx(i)),lambda(w_idx(i)),0,0);
        end
        if s
            xrange=[xrange,2*norminv(1-eps,0,s)];
        end
        xrange=max(2*xrange);
        xrange=xrange-mod(xrange,dx); % to center around 0
    end

    dx=parser.Results.dx;
    n_grid=round(parser.Results.n_grid);
    % make n_grid odd
    if mod(n_grid,2)
        n_grid=n_grid+1;
    end
    xgrid=-xrange:dx:xrange;
    xgrid=linspace(-xrange,xrange,n_grid);
    dx=2*xrange/n_grid;

    % compute non-central and normal pdf's
    % over this span
    ncpdfs=nan(nnz(w),length(xgrid));
    for i=1:nnz(w)
        pdf=gx2pdf(xgrid,w(w_idx(i)),k(w_idx(i)),lambda(w_idx(i)),0,0);
        pdf(isinf(pdf))=max(pdf(~isinf(pdf)));
        ncpdfs(i,:)=pdf;
    end
    if s
        ncpdfs=[ncpdfs;normpdf(xgrid,0,abs(s))];
    end
    f=ifft(prod(fft(ncpdfs,[],2),1));
    if any(isnan(f))||any(isinf(f))
        error('Convolution method failed. Try differentiation method.')
    end
    if ~mod(size(ncpdfs,1),2)
        f=ifftshift(f);
    end
    f=f/(sum(f)*dx);
    xgrid=xgrid+m;

    if ~strcmpi(x,'full')
        F=griddedInterpolant(xgrid,f);
        f=F(x);
    end
end