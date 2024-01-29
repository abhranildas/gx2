function f=gx2cdf_imhof_integrand(u,x,w,k,lambda,m,s)
    % define the Imhof integrand (w, k, lambda must be column vectors here)
    theta=sum(k.*atan(w*u)+(lambda.*(w*u))./(1+w.^2*u.^2),1)/2+u*(m-x)/2;
    rho=prod(((1+w.^2*u.^2).^(k/4)).*exp(((w.^2*u.^2).*lambda)./(2*(1+w.^2*u.^2))),1) .* exp(u.^2*s^2/8);
    f=sin(theta)./(u.*rho);
end