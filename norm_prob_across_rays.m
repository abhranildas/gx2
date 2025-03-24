function [p_rays,bd_pts_rays,p_tiny_sum,sym_idx]=norm_prob_across_rays(mu,v,dom,n_z,varargin)

% parse inputs
parser=inputParser;
parser.KeepUnmatched=true;
addRequired(parser,'mu',@isnumeric);
addRequired(parser,'v',@isnumeric);
addRequired(parser,'dom',@(x) isstruct(x) || isa(x,'function_handle') || ismatrix(x));
addRequired(parser,'n_z',@isnumeric);
addOptional(parser,'side','upper',@(x) strcmpi(x,'lower') || strcmpi(x,'upper') );
addParameter(parser,'output','prob',@(x) strcmpi(x,'prob') || strcmpi(x,'prob_dens') );
addParameter(parser,'dom_type','quad');
addParameter(parser,'precision','log',@(x) strcmpi(x,'basic')||strcmpi(x,'log')||strcmpi(x,'vpa'));
addParameter(parser,'fun_level',0,@isnumeric);
addParameter(parser,'fun_span',5);
addParameter(parser,'fun_resol',100);
addParameter(parser,'fun_grad',[],@(x) isa(x,'function_handle'));
addParameter(parser,'n_bd_pts',1e4);
addParameter(parser,'bd_pts',false);

parse(parser,mu,v,dom,n_z,varargin{:});

dom_type=parser.Results.dom_type;
fun_level=parser.Results.fun_level;
output=parser.Results.output;
precision=parser.Results.precision;
dim=length(mu);

n_z=n_z./vecnorm(n_z,2,1);

% ray-trace if necessary
if ~strcmpi(dom_type,'quad') || parser.Results.bd_pts
    % ray-trace the standardized domain/function
    dom_standard_raytrace=@(n) standard_ray_trace(dom,n,varargin{:},'mu',mu,'v',v);

    % initial signs and boundary distances in standardized space
    [init_sign,z]=dom_standard_raytrace(n_z);

    p_tiny_sum=0;
end

if strcmpi(dom_type,'quad')
    % standardized boundary coefficients
    quad_s=standard_quad(dom,mu,v);

    if strcmpi(precision,'basic')
        p_rays=gx2_ray_integrand(fun_level,n_z,quad_s,varargin{:});
    elseif strcmpi(precision,'log')
        [p_rays,p_tiny_sum]=gx2_ray_integrand(fun_level,n_z,quad_s,varargin{:});
    elseif strcmpi(precision,'vpa')
        [p_rays,p_tiny_sum,sym_idx]=gx2_ray_integrand(fun_level,n_z,quad_s,varargin{:});        
    end
else
    if strcmpi(output,'prob')
        % probability on rays
        if ~strcmpi(precision,'vpa')
            p_rays=cellfun(@(init_sign_ray,z_ray) prob_ray(init_sign_ray,z_ray,dim,varargin{:}),num2cell(init_sign),z);
            % if there are roots on rays with 0 prob,
            % notify to turn on vpa
            n_roots=cellfun(@(z) numel(z)>0,z);
            if nnz(n_roots&(~p_rays))
                warning("Some rays contain probabilities too small for double precision, returning 0. Set 'vpa' to true to compute these with variable precision.")
            end
        else
            p_rays=cellfun(@(init_sign_ray,z_ray) prob_ray(init_sign_ray,z_ray,dim,varargin{:}),num2cell(init_sign),z,'un',0);
        end
    elseif strcmpi(output,'prob_dens') % probability density calculations
        
        % gradient of standardized function
        fun_grad=parser.Results.fun_grad;
        standard_gradf=@(z) standard_fun_grad(z,fun_grad,mu,v);

        % probability density on rays
        p_rays=cellfun(@(n_ray,z_ray) prob_dens_ray(n_ray,z_ray,standard_gradf), num2cell(n_z,1),z);
    end

    % if p_rays is a cell but there are no symbols, convert to numeric array
    if iscell(p_rays)
        num_idx=cellfun(@isnumeric, p_rays);
        if all(num_idx)
            p_rays=cell2mat(p_rays);
        end
    end
end

if parser.Results.bd_pts
    % standard boundary points
    std_bd_pts_ray=cellfun(@(z_ray,n_ray) z_ray.*n_ray, z,num2cell(n_z,1),'un',0);

    % boundary points
    bd_pts_rays=sqrtm(v)*horzcat(std_bd_pts_ray{:})+mu;

    global bd_pts
    bd_pts=[bd_pts,bd_pts_rays];
else
    bd_pts_rays=[];
end
