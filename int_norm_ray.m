function [p,p_err,bd_pts]=int_norm_ray(mu,v,dom,varargin)
% Integrate a multinormal distribution over a specified domain, using
% the ray method.

parser = inputParser;
parser.KeepUnmatched=true;
addRequired(parser,'mu',@isnumeric);
addRequired(parser,'v',@isnumeric);
addRequired(parser,'dom');
addOptional(parser,'side','upper',@(x) strcmpi(x,'lower') || strcmpi(x,'upper') );
addParameter(parser,'dom_type','quad');
addParameter(parser,'output','prob',@(x) strcmpi(x,'prob') || strcmpi(x,'prob_dens') );
addParameter(parser,'force_mc',false,@islogical);
addParameter(parser,'fun_level',0,@isnumeric);
addParameter(parser,'fun_span',5);
addParameter(parser,'fun_resol',100);
addParameter(parser,'AbsTol',1e-10);
addParameter(parser,'RelTol',1e-2);
addParameter(parser,'vpa',false,@islogical);
addParameter(parser,'n_rays',500);
addParameter(parser,'bd_pts',false);
addParameter(parser,'gpu_batch',4e7);

parse(parser,mu,v,dom,varargin{:});
output=parser.Results.output;
fun_level=parser.Results.fun_level;
force_mc=parser.Results.force_mc;
AbsTol=parser.Results.AbsTol;
RelTol=parser.Results.RelTol;
n_rays=parser.Results.n_rays;
vpaflag=parser.Results.vpa;
dom_type=parser.Results.dom_type;
gpu_batch=parser.Results.gpu_batch;
isgpu=canUseGPU();
dim=length(mu);

if parser.Results.bd_pts
    global bd_pts
end
bd_pts=[];

if force_mc||dim>4
    % Monte-Carlo integration

    if strcmpi(dom_type,'quad')&&~vpaflag&&isgpu&&gpu_batch % batch-process on GPU
        disp('Using GPU. If this is slower, set gpu_batch to 0.');
        gpudev=gpuDevice;
        reset(gpudev)

        gpu_batch_adj=round(gpu_batch/numel(fun_level)); % adjusted GPU batch size given no. of simultaneous evaluation points
        n_batches=ceil(n_rays/gpu_batch_adj); % no. of gpu batches
        n_rays_batches=gpu_batch_adj*ones(1,n_batches); % no. of rays in each batch
        remainder=rem(n_rays,gpu_batch_adj);
        if remainder, n_rays_batches(end)=remainder; end
        p_batches=nan(numel(fun_level),n_batches);
        p2_batches=nan(numel(fun_level),n_batches); % p squared (for SEM calculation)
        for i=1:n_batches
            n_z=mvnrnd(zeros(dim,1),eye(dim),n_rays_batches(i))';
            try
                n_z=gpuArray(n_z);
                p_rays=norm_prob_across_rays(mu,v,dom,n_z,varargin{:});
                p_batches(:,i)=gather(sum(p_rays,2));
                p2_batches(:,i)=gather(sum(p_rays.^2,2));                
            catch errmsg
                warning('GPU error. If GPU is out of memory, try setting gpu_batch to a smaller value, or to 0 to not use GPU. See error below.')
                rethrow(errmsg)
            end
        end
        % mean across all batches
        p=sum(p_batches,2)/n_rays;
        % mean square across all batches
        p2=sum(p2_batches,2)/n_rays;
        % SEM of p
        p_err=sqrt((p2-p.^2)/n_rays);

    else % process all on CPU
        % uniform random rays (points on n-sphere)
        n_z=mvnrnd(zeros(dim,1),eye(dim),n_rays)';
        if ~vpaflag
            p_rays=norm_prob_across_rays(mu,v,dom,n_z,varargin{:});
            p=mean(p_rays,2);
            p_err=std(p_rays,0,2)/sqrt(n_rays);
        else % use VPA
            [p_rays,~,p_sym_sum]=norm_prob_across_rays(mu,v,dom,n_z,varargin{:});
            p_sum=sum(p_rays,2);

            % probabilities too small even with sym:
            sym_idx=cellfun(@(x) isa(x,'sym'),p_sym_sum);
            symfail_idx=logical(vpa(p_sym_sum(sym_idx))==0);
            if any(symfail_idx)
                warning('Probability on some rays too small even for variable precision. Returning 0.')
                sym_indices=find(sym_idx);
                p_sym_sum(sym_indices(symfail_idx))={0};
            end

            % include the symbolic only if it's
            % bigger than RelTol of the numeric
            log_p_sym_sum=cellfun(@(x) double(log10(x)),p_sym_sum);
            sym_shortlist_idx= sym_idx & (log_p_sym_sum > log10(p_sum) + log10(RelTol));
            if any(sym_shortlist_idx)
                warning('Probability too small for double precision. Returning as symbol, use vpa to evaluate.')
                % create cell array, combining numeric and symbolic
                % probabilitites
                p=num2cell(p_sum);
                % p(sym_shortlist_idx)=num2cell(p_sum(sym_shortlist_idx)+cell2sym(p_sym_sum(sym_shortlist_idx)));
                p(sym_shortlist_idx)=num2cell(cell2sym(p_sym_sum(sym_shortlist_idx)));
                p=cellfun(@(x) x/n_rays,p,'un',0);
                % p=(p_sum+p_sym_sum)/numel(p_rays);
                p_err=[];
            else
                p=mean(p_rays,2);
                p_err=std(p_rays,0,2)/sqrt(n_rays);
            end
        end
    end

    if strcmpi(output,'prob')
        if isnumeric(p)
            p=p/2;
        elseif iscell(p)
            p=cellfun(@(x) x/2,p,'un',0);
        end
    end

else

    % grid integration
    fun_level=parser.Results.fun_level;

    if dim==1
        p=norm_prob_across_angles(mu,v,dom,varargin{:});
    elseif dim==2
        if numel(fun_level)==1 % if integral needed only at one function level,
            % integrate with 'arrayvalued', false to evaluate
            % simultaneously across grid of theta and speed up
            p=integral(@(theta) norm_prob_across_angles(mu,v,dom,varargin{:},'theta',theta),0,pi,'AbsTol',AbsTol,'RelTol',RelTol);
        else
            % integrate with 'arrayvalued', true to simultaneously
            % integrate at multiple function levels
            p=integral(@(theta) norm_prob_across_angles(mu,v,dom,varargin{:},'theta',theta),0,pi,'arrayvalued',true,'AbsTol',AbsTol,'RelTol',RelTol);
        end
    elseif dim==3
        % use arrayfun to integrate at each fun_level separately,
        % because integral2 cannot integrate vector-valued function.
        p=arrayfun(@(f) integral2(@(theta,phi) norm_prob_across_angles(mu,v,dom,varargin{:},'theta',theta,'phi',phi,'fun_level',f),0,pi/2,0,2*pi,'AbsTol',AbsTol,'RelTol',RelTol),fun_level);
    elseif dim==4
        % use arrayfun to integrate at each fun_level separately,
        % because integral3 cannot integrate vector-valued function.
        p=arrayfun(@(f) integral3(@(theta,phi,psi) norm_prob_across_angles(mu,v,dom,varargin{:},'theta',theta,'phi',phi,'psi',psi,'fun_level',f),0,pi/2,0,2*pi,0,pi,'AbsTol',AbsTol,'RelTol',RelTol),fun_level);
    end
    p_err=[];
end