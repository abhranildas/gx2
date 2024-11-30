function [f,log10_f]=gx2pdf_tail(x,w,k,lambda,s,m,varargin)

parser=inputParser;
parser.KeepUnmatched=true;
addRequired(parser,'x',@(x) isreal(x));
addRequired(parser,'w',@(x) isreal(x) && isrow(x));
addRequired(parser,'k',@(x) isreal(x) && isrow(x));
addRequired(parser,'lambda',@(x) isreal(x) && isrow(x));
addRequired(parser,'s',@(x) isreal(x) && isscalar(x));
addRequired(parser,'m',@(x) isreal(x) && isscalar(x));
addOptional(parser,'side','upper',@(x) strcmpi(x,'lower') || strcmpi(x,'upper') );
addParameter(parser,'output','cdf',@(x) strcmpi(x,'cdf') || strcmpi(x,'pdf') );
addParameter(parser,'x_scale','linear',@(x) strcmpi(x,'linear') || strcmpi(x,'log') );

parse(parser,x,w,k,lambda,s,m,varargin{:});
side=parser.Results.side;

% first merge into unique w's
[w,~,ic]=uniquetol(w); % unique non-zero eigenvalues
k=arrayfun(@(x)sum(k(ic==x)),1:numel(w)); % merged total dof's
lambda=arrayfun(@(x)sum(lambda(ic==x)),1:numel(w)); % merged total non-centralities

if strcmpi(side,'upper') % upper tail
    [w_max,max_idx]=max(w.*(w>0));
    sgn=1;
else % lower tail
    [w_max,max_idx]=max(abs(w.*(w<0)));
    x=abs(x);
    sgn=-1;
end

k_max=k(max_idx);
lambda_max=lambda(max_idx);

w_rest=w([1:max_idx-1, max_idx+1:end]);
k_rest=k([1:max_idx-1, max_idx+1:end]);
lambda_rest=lambda([1:max_idx-1, max_idx+1:end]);

common_factor=exp(sgn*m/(2*w_max)+s^2/(8*w_max^2))*...
    prod(exp((lambda_rest.*w_rest*sgn)./(2*(w_max-sgn*w_rest)))./(1-sgn*w_rest/w_max).^(k_rest/2));

if ~lambda_max
    f=common_factor/((2*w_max)^(k_max/2)*gamma(k_max/2))*x.^(k_max/2-1).*exp(-x/(2*w_max));
else
    f=common_factor*(x/lambda_max).^(k_max/4-1/2).*exp(-x/(2*w_max)).*...
        besseli(k_max/2-1,sqrt(lambda_max*x/w_max))./...
        (2*w_max.^(k_max/4+1/2));
end

% f=common_factor*factor*x.^(k_max/2-1).*exp(-x/(2*w_max));

log10_f=(k_max/2-1)*log10(x)+log10(exp(1))*(sgn*m/(2*w_max)+s^2/(8*w_max^2)-x/(2*w_max))-...
    log10(gamma(k_max/2))-(k_max/2)*log10(2*w_max)-sum((k_rest/2).*log10(1-sgn*w_rest/w_max));