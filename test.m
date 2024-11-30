
%% debug IntClassNorm
mu_1=[4; 5]; v_1=[2 1; 1 1];
mu_2=mu_1; v_2=3*[2 -1; -1 1];

results_gx2=classify_normals([mu_1,v_1],[mu_2,v_2],'method','gx2')

%% pdf of vector function
mu=[-1; -1]; v=[1 0.5; 0.5 2];
fun=@(x,y) x.*sin(y) - y.*cos(x);

fun_grad=@(x) [sin(x(2))+x(2)*sin(x(1)); x(1)*cos(x(2))-cos(x(1))];

x=linspace(-20,20,50);
% p=norm_fun_cdf(x,mu,v,fun);
% plot(x,p)
tic
f_diff=norm_fun_pdf(x,mu,v,fun,'pdf_method','diff','dx',0.1,'AbsTol',0,'RelTol',1e-5);
time_diff=toc

tic
f_ray=norm_fun_pdf(x,mu,v,fun,'pdf_method','ray','fun_grad',fun_grad,'AbsTol',0,'RelTol',1e-5);
time_ray=toc
figure; hold on
plot(x,f_diff)
plot(x,f_ray)

%% new gx2pdf_ray
w=[-2 1 1];
k=[2 1 2];
lambda=[7 5 3];
m=5;
s=3;

x=linspace(-60,40,1e2);
tic
f=gx2pdf_ray(x,w,k,lambda,s,m,'mc_samples',1e5);
toc

plot(x,f)

% TODO near the peak of the pdf it's noisy with MC,
% and with quadrature it uses a lot of rays. Why?

%% new gx2pdf_ray, grid integral, 2d

w=[-2 1];
k=[1 1];
lambda=[7 5];
m=5;
s=0;

x=linspace(-60,40,1e2);
tic
f=gx2pdf_ray(x,w,k,lambda,s,m);
toc
% f2=gx2pdf_imhof(x,w,k,lambda,s,m); %TODO Imhof bad in this and next

figure; hold on
plot(x,f)
plot(x,f2)

%% test gx2pdf with 'method'

w=[-2 1 3];
k=[1 1 7];
lambda=[7 5 0];
m=0;
s=10;

x=linspace(-60,40,1e2);

% tic
% f_auto=gx2pdf(x,w,k,lambda,s,m);
% toc

tic
f_imhof=gx2pdf(x,w,k,lambda,s,m,'method','imhof');
toc

tic
f_ray=gx2pdf(x,w,k,lambda,s,m,'method','ray','n_rays',1e5,'vpa',false);
toc

% TODO why is ray output wiggly around -20


figure; hold on
% plot(x,f_auto)
plot(x,f_imhof,'-k')
plot(x,f_ray,'-b')
% plot(x,f_conv)

%% gx2inv with log prob

w=[-2 1 3];
k=[1 1 7];
lambda=[7 5 0];
m=0;
s=10;

p=-400;

x=gx2inv([p p],w,k,lambda,s,m,'method','ray','vpa',true,'mc_samples',5e4);



%% test tail pdf and cdf approximation

w=[-2 -3 4 1];
k=[4 2 6 1];
lambda=[10 20 5 0];
s=20;
m=0;

x=linspace(0,4e3,20);

% upper tail

% pdf
f_ray=gx2pdf(x,w,k,lambda,s,m,'method','ray','n_rays',1e7,'vpa',false);
f_imhof=gx2pdf(x,w,k,lambda,s,m,'method','imhof','vpa',false,'abstol',1e-100,'reltol',1e-100);
f_pearson=gx2pdf(x,w,k,lambda,s,m,'method','pearson');
f_tail=gx2pdf(x,w,k,lambda,s,m,'upper','method','das');

figure; hold on
plot(x,log10(f_imhof),'-ok')
plot(x,log10(f_ray),'-ob')
plot(x,log10(f_pearson),'-og')
plot(x,log10(f_tail),'-or')

% cdf
p_ray=gx2cdf(x,w,k,lambda,s,m,'upper','method','ray','n_rays',1e7,'vpa',false);
p_imhof=gx2cdf(x,w,k,lambda,s,m,'upper','method','imhof','vpa',false,'abstol',1e-100,'reltol',1e-100);
p_pearson=gx2cdf(x,w,k,lambda,s,m,'upper','method','pearson');
p_tail=gx2cdf(x,w,k,lambda,s,m,'upper','method','das');

figure; hold on
plot(x,log10(p_imhof),'-ok')
plot(x,log10(p_ray),'-ob')
plot(x,log10(p_pearson),'-og')
plot(x,log10(p_tail),'-or')

% lower tail

x=-linspace(0,4e3,20);

% pdf
f_ray=gx2pdf(x,w,k,lambda,s,m,'method','ray','n_rays',1e7,'vpa',false);
f_imhof=gx2pdf(x,w,k,lambda,s,m,'method','imhof','vpa',false,'abstol',1e-200,'reltol',1e-200);
f_pearson=gx2pdf(x,w,k,lambda,s,m,'method','pearson');
f_tail=gx2pdf(x,w,k,lambda,s,m,'method','das');

figure; hold on
plot(x,log10(f_imhof),'-ok')
plot(x,log10(f_ray),'-b')
plot(x,log10(f_pearson),'-g')
plot(x,log10(f_tail),'-r')

% cdf
p_ray=gx2cdf(x,w,k,lambda,s,m,'method','ray','n_rays',1e7,'vpa',false);
p_imhof=gx2cdf(x,w,k,lambda,s,m,'method','imhof','vpa',false,'abstol',1e-200,'reltol',1e-200);
p_pearson=gx2cdf(x,w,k,lambda,s,m,'method','pearson');
p_tail=gx2cdf(x,w,k,lambda,s,m,'method','das');

figure; hold on
plot(x,log10(p_imhof),'-ok')
plot(x,log10(p_ray),'-b')
plot(x,log10(p_pearson),'-g')
plot(x,log10(p_tail),'-r')

%% log sum exp approximation for ray vpa

w=[-2 -3 4 1];
k=[4 2 6 1];
lambda=[10 20 5 0];
s=20;
m=0;

x=linspace(0,1e4,20);

% cdf
[p_ray_log,p_err]=gx2cdf(x,w,k,lambda,s,m,'upper','method','ray','n_rays',1e7);


figure; hold on
plot(x(p_ray_log>0),log10(p_ray_log(p_ray_log>0)),'-ob','MarkerFaceColor','w')
plot(x(p_ray_log<0),p_ray_log(p_ray_log<0),'-ob','MarkerFaceColor','b')
% plot(x,double(log10(cell2sym(p_ray))),'-ob')

%% log sum exp approximation for ray vpa

w=[2 3 4 1];
k=[4 2 6 1];
lambda=[0 0 0 0];
s=0;
m=0;

x=logspace(-70,-40,20);

% cdf
tic
p_ray_log=gx2cdf(x,w,k,lambda,s,m,'method','ray','n_rays',1e6);
toc
% tic
% p_ray_vpa=gx2cdf(x,w,k,lambda,s,m,'method','ray','precision','vpa','n_rays',10);
% toc

figure; hold on
plot(x(p_ray_log>0),log10(p_ray_log(p_ray_log>0)),'-ob','MarkerFaceColor','w')
plot(x(p_ray_log<0),p_ray_log(p_ray_log<0),'-ob','MarkerFaceColor','b')
set(gca,'xscale','log')
% plot(x,double(log10(cell2sym(p_ray_vpa))),'-or')

