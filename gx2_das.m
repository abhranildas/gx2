function [p,errflag]=gx2_das(x,w,k,lambda,s,m,varargin)

parser=inputParser;
parser.KeepUnmatched=true;
addRequired(parser,'x',@(x) isreal(x));
addRequired(parser,'w',@(x) isreal(x) && isrow(x));
addRequired(parser,'k',@(x) isreal(x) && isrow(x));
addRequired(parser,'lambda',@(x) isreal(x) && isrow(x));
addRequired(parser,'s',@(x) isreal(x) && isscalar(x));
addRequired(parser,'m',@(x) isreal(x) && isscalar(x));
addOptional(parser,'side','lower',@(x) strcmpi(x,'lower') || strcmpi(x,'upper') );
addParameter(parser,'output','cdf',@(x) strcmpi(x,'cdf') || strcmpi(x,'pdf') );
addParameter(parser,'scale','linear',@(x) strcmpi(x,'linear') || strcmpi(x,'log') );

parse(parser,x,w,k,lambda,s,m,varargin{:});
side=parser.Results.side;

% first merge into unique w's
[w,~,ic]=uniquetol(w); % unique non-zero eigenvalues
k=arrayfun(@(x)sum(k(ic==x)),1:numel(w)); % merged total dof's
lambda=arrayfun(@(x)sum(lambda(ic==x)),1:numel(w)); % merged total non-centralities

if strcmpi(side,'upper') % upper tail
    [w_max,max_idx]=max(w.*(w>0));
else % lower tail
    [w_max,max_idx]=min(w.*(w<0));
    % x=abs(x);
    % sgn=-1;
end
sgn=1;

k_max=k(max_idx);
lambda_max=lambda(max_idx);

w_rest=w([1:max_idx-1, max_idx+1:end]);
k_rest=k([1:max_idx-1, max_idx+1:end]);
lambda_rest=lambda([1:max_idx-1, max_idx+1:end]);

a=exp(sgn*m/(2*w_max)+s^2/(8*w_max^2))*...
    prod(exp((lambda_rest.*w_rest*sgn)./(2*(w_max-sgn*w_rest)))./(1-sgn*w_rest/w_max).^(k_rest/2));

if strcmpi(parser.Results.output,'pdf')
    if strcmpi(parser.Results.scale,'linear')
        p=a/abs(w_max)*ncx2pdf(x/w_max,k_max,lambda_max);
    % elseif strcmpi(parser.Results.scale,'linear')
        % p=(k_max/2-1)*log10(x)+log10(exp(1))*(sgn*m/(2*w_max)+s^2/(8*w_max^2)-x/(2*w_max))-...
    % log10(gamma(k_max/2))-(k_max/2)*log10(2*w_max)-sum((k_rest/2).*log10(1-sgn*w_rest/w_max));
    end
elseif strcmpi(parser.Results.output,'cdf')
    marcum=marcumq(sqrt(lambda_max),sqrt(x/w_max),k_max/2);
    p=a*marcum;
end

errflag=p<0;

if any(errflag)
    warning('Tail approximation output(s), i.e. Marcum Q-function values, are too small to correctly compute, so clipping to 0. Check the error flag output.')
    p=max(p,0);
end

log10_f=(k_max/2-1)*log10(x)+log10(exp(1))*(sgn*m/(2*w_max)+s^2/(8*w_max^2)-x/(2*w_max))-...
    log10(gamma(k_max/2))-(k_max/2)*log10(2*w_max)-sum((k_rest/2).*log10(1-sgn*w_rest/w_max));