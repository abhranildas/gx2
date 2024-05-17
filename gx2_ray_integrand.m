function [p_rays,p_sym_sum,sym_idx]=gx2_ray_integrand(x,n_z,quad,varargin)
% return the differential probability or probability density on each ray
% that is integrated across rays

% parse inputs
parser=inputParser;
parser.KeepUnmatched=true;
addRequired(parser,'x',@isnumeric); % points at which to find the cdf/pdf
addRequired(parser,'n_z',@isnumeric);
addRequired(parser,'quad',@isstruct);
addOptional(parser,'side','lower',@(x) strcmpi(x,'lower') || strcmpi(x,'upper') );
addParameter(parser,'output','prob'); % probability or probability density
addParameter(parser,'vpa',false,@islogical);

parse(parser,x,n_z,quad,varargin{:});

output=parser.Results.output;
side=parser.Results.side;
vpaflag=parser.Results.vpa;

dim=numel(quad.q1);

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
    [f_big,f_small]=Phibar_ray_split(z,dim);
    p_rays_big=init_sign_rays+1+init_sign_rays.*(f_big(:,:,2)-f_big(:,:,1));
    p_rays_small=init_sign_rays.*(f_small(:,:,2)-f_small(:,:,1));
    if strcmpi(side,'upper')
        p_rays=p_rays_big+p_rays_small;
    elseif strcmpi(side,'lower')
        p_rays=2-p_rays_big-p_rays_small;
    end
elseif strcmpi(output,'prob_dens')

    % phi_ray at each root
    sum_phi=sum(phi_ray(z,dim),3,'omitnan');

    % quadratic slope at each root
    % quad_slope=nan(numel(x),size(n_z,2));
    % quad_slope(quad_root_exists)=delta(quad_root_exists);
    % quad_slope(linear_root_exists)=abs(q1(linear_root_exists));
    quad_slope=nan(size(delta2));
    quad_slope(root_exists)=sqrt(delta2(root_exists)); % slope is the same formula at quad and linear roots

    % divide phi_ray by quad slope
    p_rays=sum_phi./quad_slope;

    p_rays(isnan(p_rays))=0;
end

p_sym_sum=num2cell(zeros(numel(x),1));
sym_idx=[];
% cases where root exists but computed prob. is 0,
tiny_probs=root_exists&(~p_rays);
if nnz(tiny_probs)
    if ~vpaflag
        warning("%.1f%% of rays contain probabilities smaller than realmin=1e-308, returning 0. Set 'vpa' to true to compute these with variable precision.",100*mean(tiny_probs,'all'))
    else
        sym_idx=any(tiny_probs,2); % levels at which at least one ray is sym.
        % group roots into a cell array so we can use cellfun:
        z_cell=permute(num2cell(permute(z,[1 3 2]),2),[1 3 2]);
        % remove nan roots:
        z_cell=cellfun(@(z_each) z_each(~isnan(z_each)),z_cell,'un',0);
        if strcmpi(output,'prob')
            p_sym_sum(sym_idx)=sym2cell(sum(cellfun(@(init_sign_ray,z_ray) prob_ray_sym(init_sign_ray,z_ray,dim,side),num2cell(init_sign_rays(sym_idx,:)),z_cell(sym_idx,:)),2));
        elseif strcmpi(output,'prob_dens')
            % p_sym_sum(sym_idx)=sym2cell(arrayfun(@(iLevel) ...
            %     sum(sum(phi_ray(sym(z(iLevel,tiny_probs(iLevel,:),:)),dim),3)./quad_slope(iLevel,tiny_probs(iLevel,:))),...
            %     find(sym_idx)));
            p_sym=cellfun(@(z_ray,quad_slope_ray) sum(phi_ray(sym(z_ray),dim))/quad_slope_ray,z_cell(sym_idx,:),num2cell(quad_slope(sym_idx,:)));
            p_sym_sum(sym_idx)=sym2cell(sum(p_sym,2,'omitnan'));
        end
    end
end