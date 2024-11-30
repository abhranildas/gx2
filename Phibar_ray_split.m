function [Phibar_big,Phibar_small]=Phibar_ray_split(z,dim)
    % returns the complementary cdf of the standard multinormal on a ray, split into a big chunk, which is
    % 0, +/-2 or +/-1, and a small part, which is a tail CDF. This is to
    % prevent the small part vanishing when the two are added together due to machine imprecision.

    z_c=sqrt(chi2inv(0.5,dim)); % if z^2> z_c^2, we'll use chi2cdf('upper').
    Phibar_big=(z<=-z_c)+(z<z_c); % 2 when z<=-z_c, 1 when -z_c<z<z_c, 0 when z>=z_c

    Phibar_small=gammainc(z.^2/2,dim/2); %chi2cdf(z.^2,dim);
    Phibar_small_upper=gammainc(z.^2/2,dim/2,'upper'); %chi2cdf(z.^2,dim,'upper');
    Phibar_small(abs(z)>=z_c)=Phibar_small_upper(abs(z)>=z_c);

    Phibar_small=Phibar_small.*(1-2*((z<=-z_c)|((z>0)&(z<z_c)))); % if z<=-z_c or 0<z<z_c, invert sign

    % if any z is nan (non-existent root),
    % prob. contribution is 0.
    Phibar_big(isnan(Phibar_big))=0;
    Phibar_small(isnan(Phibar_small))=0;

end