function p=Phibar_sym(z,dim)
% complementary cdf on ray using symbol-friendly igamma
% first use igamma to construct regularized gamma, using which
% construct chi distribution cdf, using which construct phibar.

% chi CDF is regularized gamma, built from symbol-friendly igamma
chiCDF=1-igamma(dim/2,z.^2/2)/gamma(dim/2);

% Phibar from chicdf
p=(1-sign(z).*chiCDF)/2;
end