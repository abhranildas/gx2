function f=phi_ray(z,dim)
    % standard multinormal density function along a ray
    % f=abs(z).*chi2pdf(z.^2,dim);
    z2=z.^2;
    f=abs(z).*(z2.^(dim/2-1)).*exp(-z2/2)/(2^(dim/2)*gamma(dim/2));
end