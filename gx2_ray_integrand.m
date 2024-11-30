function [p_rays,p_tiny_sum,sym_idx]=gx2_ray_integrand(x,n_z,quad,varargin)
% return the differential probability or probability density on each ray
% that is integrated across rays

% parse inputs
parser=inputParser;
parser.KeepUnmatched=true;
addRequired(parser,'x',@isnumeric); % points at which to find the cdf/pdf
addRequired(parser,'n_z',@isnumeric);
addRequired(parser,'quad',@isstruct);
addOptional(parser,'side','upper',@(x) strcmpi(x,'lower') || strcmpi(x,'upper') );
addParameter(parser,'output','prob'); % probability or probability density
addParameter(parser,'precision','log',@(x) strcmpi(x,'basic')||strcmpi(x,'log')||strcmpi(x,'vpa'));

parse(parser,x,n_z,quad,varargin{:});

output=parser.Results.output;
side=parser.Results.side;
precision=parser.Results.precision;

dim=numel(quad.q1);
n_levels=numel(x);

% find the quadratic coefficients across all rays
q2=dot(n_z,quad.q2*n_z,1);
q1=quad.q1'*n_z;
q0=quad.q0;

% discriminant of the quadratic across rays and levels
delta2=q1.^2-4*q2.*(q0-x); % delta^2
root_exists=delta2>0; % levels where linear or quadratic roots exist
quad_root_exists=root_exists & q2; % levels where quadratic roots exist
delta=nan(size(delta2));
delta(quad_root_exists)=sqrt(delta2(quad_root_exists)); % populate only with quad delta for now

% linear_root_exists=repmat(~q2 & q1, [numel(x) 1]); % levels where linear roots exist

% sorted roots across rays
z=(-q1+cat(3,-1,1).*delta)./(2*abs(q2)); % quadratic roots where q2 ~= 0
if nnz(~q2)
    z(:,~q2,1)=-(q0-x)./q1(~q2); % linear roots where q2=0
end

if strcmpi(output,'prob')
    init_sign_rays=sign(4*sign(q2)-2*sign(q1)+sign(q0-x));
    [Phibar_big,Phibar_small]=Phibar_ray_split(z,dim);
    p_rays_big=init_sign_rays+1+init_sign_rays.*(Phibar_big(:,:,2)-Phibar_big(:,:,1));
    p_rays_small=init_sign_rays.*(Phibar_small(:,:,2)-Phibar_small(:,:,1));
    if strcmpi(side,'upper')
        p_rays=p_rays_big+p_rays_small;
    elseif strcmpi(side,'lower')
        p_rays=2-p_rays_big-p_rays_small;
    end
elseif strcmpi(output,'prob_dens')
    % phi_ray at each root
    sum_phi=sum(phi_ray(z,dim),3,'omitnan');

    % quadratic slope at each root
    quad_slope=nan(size(delta2));
    quad_slope(root_exists)=sqrt(delta2(root_exists)); % slope is the same formula at quad and linear roots

    % divide phi_ray by quad slope
    p_rays=sum_phi./quad_slope;

    p_rays(isnan(p_rays))=0;
end

% cases where root exists but computed prob. is 0,
tiny_probs=root_exists&(~p_rays);

if strcmpi(precision,'basic')
    if nnz(tiny_probs)
        warning("%.1f%% of rays contain probabilities less than realmin=1e-308, returning 0. Set 'precision' to 'log' or 'vpa' to compute these.",100*mean(tiny_probs,'all'))
    end
elseif strcmpi(precision,'log')
    % log sum exp of tiny probabilities
    if strcmpi(output,'prob')
        % roots where Phibar_small was tiny, or root pairs where they weren't tiny
        % individually but cancelled
        % tiny_probs=repmat(p_rays==0,[1 1 2])|(Phibar_small==0);
        tiny_probs=repmat(tiny_probs,[1 1 2]);
        z_tiny=nan(size(z));
        z_tiny(tiny_probs)=z(tiny_probs);

        z_tiny_signs=init_sign_rays.*sign(z_tiny).*cat(3,-1,1); % sign of contribution from each root
        % z_tiny_signs=nan(size(z));
        % z_tiny_signs(tiny_probs&~isnan(z_tiny))=all_signs(tiny_probs&~isnan(z_tiny));

        % log of tiny Phibar smalls:
        log_Phibar_small=(dim-2)*log10(abs(z_tiny))-z_tiny.^2/(2*log(10))-log10(gamma(dim/2)*2^(dim/2-1));

        % log sum exp of tiny Phibar smalls with positive and negative signs:
        log_Phibar_plus=nan(n_levels,1);
        log_Phibar_minus=nan(n_levels,1);
        for level=1:n_levels
            log_Phibar_small_thislevel=log_Phibar_small(level,:,:);
            log_Phibar_plus(level)=log_sum_exp(log_Phibar_small_thislevel(z_tiny_signs(level,:,:)==1));
            log_Phibar_minus(level)=log_sum_exp(log_Phibar_small_thislevel(z_tiny_signs(level,:,:)==-1));
        end

        % now subtract minus from plus:
        % 1. sign:
        p_tiny_sign=2*((log_Phibar_plus>log_Phibar_minus)-.5);

        % 2. magnitude:

        % first find max and min of each pair element-wise
        max_log=max(log_Phibar_plus,log_Phibar_minus);
        min_log=min(log_Phibar_plus,log_Phibar_minus);

        % then compute log10(|a - b|) using the formula
        log_Phibar_abs=max_log + log10(abs(1 - 10.^(min_log - max_log)));
        log_Phibar_abs(min_log==max_log)=-inf;
        p_tiny_sum=p_tiny_sign.*log_Phibar_abs; % combine sign with the log
        
    end

elseif strcmpi(precision,'vpa')
    p_tiny_sign=[];
    p_tiny_sum=num2cell(zeros(n_levels,1));
    sym_idx=[];
    if nnz(tiny_probs)
        sym_idx=any(tiny_probs,2); % levels at which at least one ray is sym.
        % group roots into a cell array so we can use cellfun:
        z_cell=permute(num2cell(permute(z,[1 3 2]),2),[1 3 2]);
        % remove nan roots:
        z_cell=cellfun(@(z_each) z_each(~isnan(z_each)),z_cell,'un',0);
        if strcmpi(output,'prob')
            p_tiny_sum(sym_idx)=sym2cell(sum(cellfun(@(init_sign_ray,z_ray) prob_ray_sym(init_sign_ray,z_ray,dim,side),num2cell(init_sign_rays(sym_idx,:)),z_cell(sym_idx,:)),2));
        elseif strcmpi(output,'prob_dens')
            p_sym=cellfun(@(z_ray,quad_slope_ray) sum(phi_ray(sym(z_ray),dim))/quad_slope_ray,z_cell(sym_idx,:),num2cell(quad_slope(sym_idx,:)));
            p_tiny_sum(sym_idx)=sym2cell(sum(p_sym,2,'omitnan'));
        end
    end
end