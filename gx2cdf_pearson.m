function [p,flag]=gx2cdf_pearson(x,w,k,lambda,m,varargin)
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