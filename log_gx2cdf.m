function l=log_gx2cdf(x,w,k,lambda,s,m,varargin)
p=gx2cdf(x,w,k,lambda,s,m,varargin{:});
if isnumeric(p)
    if p<=0
        l=p;
    elseif p>0
        l=log10(p);
    end
elseif iscell(p)
    l=double(log10(cell2sym(p)));
end
% end
end