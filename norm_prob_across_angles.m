function [p_angles,bd_pts_angles]=norm_prob_across_angles(mu,v,dom,varargin)

    % parse inputs
    parser=inputParser;
    parser.KeepUnmatched=true;
    addRequired(parser,'mu',@isnumeric);
    addRequired(parser,'v',@isnumeric);
    addRequired(parser,'dom',@(x) isstruct(x) || isa(x,'function_handle') || ismatrix(x));
    addOptional(parser,'side','upper',@(x) strcmpi(x,'lower') || strcmpi(x,'upper') );
    addParameter(parser,'dom_type','quad');
    addParameter(parser,'output','prob'); % probability or probability density
    addParameter(parser,'fun_span',5);
    addParameter(parser,'fun_resol',100);
    addParameter(parser,'n_bd_pts',500);
    addParameter(parser,'bd_pts',false);
    addParameter(parser,'theta',nan,@isnumeric);
    addParameter(parser,'phi',nan,@isnumeric);
    addParameter(parser,'psi',nan,@isnumeric);

    parse(parser,mu,v,dom,varargin{:});

    n_bd_pts=parser.Results.n_bd_pts;
    theta=parser.Results.theta;
    phi=parser.Results.phi;
    psi=parser.Results.psi;
    output=parser.Results.output;

    dim=length(mu);
    if dim>1 && any(strcmpi(parser.UsingDefaults,'theta')) && any(strcmpi(parser.UsingDefaults,'phi'))
        points=fibonacci_sphere(n_bd_pts);
        [az,el]=cart2sph(points(1,:),points(2,:),points(3,:));
        theta=pi/2-el;
        phi=az;
        if dim==2
            theta=linspace(0,pi,n_bd_pts);
        end
    end

    if dim==1
        n_z=1;
    elseif dim==2
        [x,y]=pol2cart(theta,1);
        n_z=[x;y];
    elseif dim==3
        az=phi;
        el=pi/2-theta;
        [x,y,z]=sph2cart(az,el,1);
        n_z=[x(:)';y(:)';z(:)'];
    elseif dim==4
        x=cos(psi);
        y=sin(psi).*cos(theta);
        z=sin(psi).*sin(theta).*cos(phi);
        w=sin(psi).*sin(theta).*sin(phi);
        n_z=[x(:)';y(:)';z(:)';w(:)'];
    end

    [p_angles,bd_pts_angles]=norm_prob_across_rays(mu,v,dom,n_z,varargin{:});

    % factors: first multiply by 2/total angle
    if dim==2
        p_angles=p_angles/pi;
    elseif dim==3
        p_angles=reshape(p_angles,size(theta));
        p_angles=p_angles.*sin(theta)/(2*pi);
    elseif dim==4
        p_angles=reshape(p_angles,size(theta));
        p_angles=p_angles.*sin(theta).*sin(psi).^2/(pi^2);
    end

    if strcmpi(output,'prob')
        % for prob. density, the factor is 2/total angle,
        % but for prob, it's 1/total angle.
        % so extra factor of 1/2.
        p_angles=p_angles/2;
    end