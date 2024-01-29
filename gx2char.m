function phi=gx2char(t,w,k,lambda,m,s)

    parser = inputParser;
    addRequired(parser,'w',@(x) isreal(x) && isrow(x));
    addRequired(parser,'k',@(x) isreal(x) && isrow(x));
    addRequired(parser,'lambda',@(x) isreal(x) && isrow(x));
    addRequired(parser,'m',@(x) isreal(x) && isscalar(x));
    addRequired(parser,'s',@(x) isreal(x) && isscalar(x));
    parse(parser,w,k,lambda,m,s);

    t_flat=t(:); % flatten input array x into a column vector of levels y

    phi=exp(1i*m*t_flat+1i*t_flat.*sum((w.*lambda)./(1-2i*t_flat.*w),2)-s^2*t_flat.^2/2)./...
        prod((1-2i*w.*t_flat).^(k/2),2);

    % reshape output to shape of input t
    phi=reshape(phi,size(t));

