function l=log_gx2cdf(x,w,k,lambda,s,m,varargin)
p=gx2cdf(x,w,k,lambda,s,m,varargin{:});

if p<=0
    l=p;
elseif p>0
    l=log10(p);
end
