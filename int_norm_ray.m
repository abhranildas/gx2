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
addParameter(parser,'precision','log',@(x) strcmpi(x,'basic')||strcmpi(x,'log')||strcmpi(x,'vpa'));
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
precision=parser.Results.precision;
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

    if strcmpi(dom_type,'quad') && ~strcmpi(precision,'vpa') && isgpu && gpu_batch % batch-process on GPU
        disp('Using GPU. To use only CPU, set gpu_batch to 0.');
        gpudev=gpuDevice;
        reset(gpudev)

        gpu_batch_adj=round(gpu_batch/numel(fun_level)); % adjusted GPU batch size given no. of simultaneous evaluation points
        n_batches=ceil(n_rays/gpu_batch_adj); % no. of gpu batches
        n_rays_batches=gpu_batch_adj*ones(1,n_batches); % no. of rays in each batch
        remainder=rem(n_rays,gpu_batch_adj);
        if remainder, n_rays_batches(end)=remainder; end
        p_batches=nan(numel(fun_level),n_batches);
        p2_batches=nan(numel(fun_level),n_batches); % p squared (for SEM calculation)
        p_tiny_sum_batches=nan(numel(fun_level),n_batches); % for log precision
        for i=1:n_batches
            n_z=mvnrnd(zeros(dim,1),eye(dim),n_rays_batches(i))';
            try
                n_z=gpuArray(n_z);
                if strcmpi(precision,'basic')
                    p_rays=norm_prob_across_rays(mu,v,dom,n_z,varargin{:});
                elseif strcmpi(precision,'log')
                    [p_rays,~,p_tiny_sum_batches(:,i)]=norm_prob_across_rays(mu,v,dom,n_z,varargin{:});
                end
                p_batches(:,i)=gather(sum(p_rays,2));
                p2_batches(:,i)=gather(sum(p_rays.^2,2));
            catch errmsg
                warning('GPU error. If GPU is out of memory, try setting gpu_batch to a smaller value, or to 0 to not use GPU. See error below.')
                rethrow(errmsg)
            end
        end
        % mean across all batches
        p_sum=sum(p_batches,2);
        p=p_sum/n_rays;
        % mean square across all batches
        p2=sum(p2_batches,2)/n_rays;
        % SEM of p
        p_err=sqrt((p2-p.^2)/n_rays);

        if strcmpi(precision,'log')
            p_tiny_sum=signed_log_sum_exp(p_tiny_sum_batches,2);
        end

    else % process all on CPU
        % uniform random rays (points on n-sphere)
        n_z=mvnrnd(zeros(dim,1),eye(dim),n_rays)';
        if strcmpi(precision,'basic')
            p_rays=norm_prob_across_rays(mu,v,dom,n_z,varargin{:});
            p=mean(p_rays,2);
            p_err=std(p_rays,0,2)/sqrt(n_rays);
        elseif strcmpi(precision,'log')
            [p_rays,~,p_tiny_sum]=norm_prob_across_rays(mu,v,dom,n_z,varargin{:});
            p_sum=sum(p_rays,2);

            p_err=std(p_rays,0,2)/sqrt(n_rays);

        elseif strcmpi(precision,'vpa')
            [p_rays,~,p_tiny_sum,sym_idx]=norm_prob_across_rays(mu,v,dom,n_z,varargin{:});
            p_sum=sum(p_rays,2);

            % probabilities too small even with sym:
            symfail_idx=logical(vpa(p_tiny_sum(sym_idx))==0);
            if any(symfail_idx)
                warning('Probability on some rays too small even for variable precision. Returning 0.')
                sym_indices=find(sym_idx);
                p_tiny_sum(sym_indices(symfail_idx))={0};
            end

            % include the symbolic only if it's
            % bigger than RelTol of the numeric
            log_p_sym_sum=cellfun(@(x) double(log10(x)),p_tiny_sum);
            sym_shortlist_idx= sym_idx & (log_p_sym_sum > log10(p_sum) + log10(RelTol));
            if any(sym_shortlist_idx)
                warning('Probability too small for double precision. Returning as symbol, use vpa to evaluate.')
                % create cell array, combining numeric and symbolic
                % probabilitites
                p_sum=num2cell(p_sum);
                p_sum(sym_shortlist_idx)=num2cell(cell2sym(p_tiny_sum(sym_shortlist_idx)));
                p=cellfun(@(x) x/n_rays,p_sum,'un',0);
                p_err=[];
            else
                p=mean(p_rays,2);
                p_err=std(p_rays,0,2)/sqrt(n_rays);
            end
        end
    end

    if strcmpi(precision,'log')
        % merge the log with the main, correctly interpreting the sign
        % of the log
        p_tiny_sign=-sign(p_tiny_sum);
        p_sum=p_sum+p_tiny_sign.*(10.^(-abs(p_tiny_sum)));
        p=p_sum/n_rays;
        % p_err=p_tiny_sum;
        if any(p_tiny_sign(~p)==-1)
            error('p_tiny has wrong sign')
        end

        % tiny_idx=(~p)&(p_tiny_sum>-inf); % places where p=0, but p_tiny has some value
        % p(tiny_idx)=p_tiny_sum(tiny_idx)-log10(n_rays);
        p_tiny_sum=p_tiny_sum-log10(n_rays); % divide by the number of rays
        p_tiny_sum(p_tiny_sum==-inf)=0; % set exactly 0 to 0
        p(~p)=p_tiny_sum(~p);

        if any(p<0)
            warning('Some output values are small for double precision. Returning their log10 values, which are negative.')
        end
    end

    % divide probs by 2
    if strcmpi(output,'prob')
        if isnumeric(p)
            p(p>0)=p(p>0)/2;
            p(p<0)=p(p<0)-log10(2);
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