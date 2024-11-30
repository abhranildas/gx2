function p_ray=prob_ray_sym(init_sign_ray,z_ray,dim,side)

    % function to compute probability slice along each ray symbolically, for
    % tiny probabilities

    sgns=sym(init_sign_ray);
    dim_s=sym(dim);

    Phibars=arrayfun(@(z) Phibar_sym(z,dim_s),z_ray);

    p_ray=sgns+1+2*sgns*sum((-1).^(1:length(z_ray)).*Phibars);

    if exist('side','var') && strcmpi(side,'lower')
        p_ray=2-p_ray;
    end

end