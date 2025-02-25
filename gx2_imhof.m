function [p,errflag]=gx2_imhof(x,w,k,lambda,s,m,varargin)

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
addParameter(parser,'precision','basic',@(x) strcmpi(x,'basic')||strcmpi(x,'vpa'));
addParameter(parser,'AbsTol',1e-10,@(x) isreal(x) && isscalar(x) && (x>=0));
addParameter(parser,'RelTol',1e-2,@(x) isreal(x) && isscalar(x) && (x>=0));

parse(parser,x,w,k,lambda,s,m,varargin{:});
output=parser.Results.output;
side=parser.Results.side;
AbsTol=parser.Results.AbsTol;
RelTol=parser.Results.RelTol;

% compute the integral
if strcmpi(parser.Results.precision,'basic')
    imhof_integral=arrayfun(@(x) integral(@(u) gx2_imhof_integrand(u,x,w',k',lambda',s,m,output),0,inf,'AbsTol',AbsTol,'RelTol',RelTol),x);
elseif strcmpi(parser.Results.precision,'vpa')
    syms u
    imhof_integral=arrayfun(@(x) vpaintegral(@(u) gx2_imhof_integrand(u,x,w',k',lambda',s,m,output),...
        u,0,inf,'AbsTol',AbsTol,'RelTol',RelTol),x);
end

if strcmpi(output,'cdf')
    if strcmpi(side,'lower')
        p=0.5-imhof_integral/pi;
    elseif strcmpi(side,'upper')
        p=0.5+imhof_integral/pi;
    end
    if isa(p,'sym')
        p=double(vpa(p));
    end
    % if p_double>realmin
        % p=p_double;
    % end
    % if isa(p, 'double')
        errflag = p<0 | p>1;
        p=min(p,1);
    % end
elseif strcmpi(output,'pdf')
    p=imhof_integral/(2*pi);
    if isa(p,'sym')
        p=double(vpa(p));
    end
    % p_double=double(vpa(p));
    % if p_double>realmin
        % p=p_double;
    % end
    % if isa(p, 'double')
        errflag = p<0;
    % end
end

if any(errflag)
    warning('Imhof method output(s) too close to limit to compute exactly, so clipping. Check the flag output, and try stricter tolerances.')
    p=max(p,0);
end

end