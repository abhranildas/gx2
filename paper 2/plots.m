%% define colours
orange=[255 123 49]/255;
sky=[0 179 255]/255;
cobalt=[0 116 255]/255;

%% quadratic picture corresponding to gx2cdf

w=[1 -1];
k=[1 1];
lambda=[2 4];
m=0;
s=0;

x_full=[linspace(-20,-.5,200) linspace(-.5,.5,100) linspace(.5,20,50)];
x=linspace(-15,15,7);

colors=spring(length(x));

% plot pdf
f_full=gx2pdf(x_full,w,k,lambda,m,s,'method','conv');
figure; hold on

int_idx=x_full<x(3);
a=area(x_full(int_idx),f_full(int_idx));
a(1).FaceColor=.1*blue+.9*[1 1 1];%[.9,1,.96];
for i=1:length(x)
    xline(x(i),'color',colors(i,:),'linewidth',.75)
end
plot(x_full,f_full,'color',blue,'linewidth',1)
xlim([-20 20])
set(gca,'xtick',[-20 0 20],'ytick',[],'fontsize',13)
xlabel("$\tilde{\chi}$",'interpreter','latex')
ylabel 'pdf'

% plot normal probability view

quad=gx2_to_norm_quad_params(w,k,lambda,m,s);
% flip quad for lower side
quad=structfun(@uminus,quad,'un',0);
mu=[0;0]; v=eye(2);

figure; hold on
plot_normal(mu,v,1,blue)

quad_this=quad;
quad_this.q0=quad_this.q0+x(3);
plot_boundary(quad_this,2,'plot_type','fill','fill_colors',blue);

p=nan(size(x));
for i=1:length(x)
    quad_this=quad;
    quad_this.q0=quad_this.q0+x(i);
    p(i)=integrate_normal(mu,v,quad_this,'method','ray','plotmode',0);
    plot_boundary(quad_this,2,'plot_type','line','line_color',colors(i,:));
end
axis([-5 5 -5 5])
xlabel("$z_1$",'interpreter','latex')
ylabel("$z_2$",'interpreter','latex')

set(gca,'xtick',[],'ytick',[],'fontsize',13)

% plot cdf
p_full=gx2cdf(x_full,w,k,lambda,m,s);
figure; hold on
plot(x_full,p_full,'color',blue,'linewidth',.75)
for i=1:length(x)
    plot(x(i),p(i),'o','markerfacecolor',colors(i,:),'markeredgecolor','k','markersize',4)
end
set(gca,'xtick',[],'ytick',[0 1],'fontsize',13)
ylabel 'cdf'

%% compare gx2 methods with ray method
% just an ncx2

w=1;
k=4;
lambda=5;
m=0;
s=0;

% find support
med=gx2inv(0.5,w,k,lambda,m,s);
[~,v]=gx2stat(w,k,lambda,m,s);

% bottom end
x_bot=logspace(-170,log10(.08),10);
x_bot=x_bot(1:end-1);

% middle
x_midbot=linspace(.08,med,10);
x_midtop=linspace(med,25,13);

% top end
x_top=[logspace(log10(25),4,20) logspace(7,9,5)];

% ncx2 method
tic
p_bot_ncx2=gx2cdf([x_bot x_midbot],w,k,lambda,m,s);
p_top_ncx2=gx2cdf([x_midtop x_top],w,k,lambda,m,s,'upper');
time_ncx2=toc

% ray method
% tic
[p_bot_ray,t_bot_ray,vpaflag_bot_ray]=gx2cdf_ray([x_bot x_midbot],w,k,lambda,m,s,'mc_samples',2e3);
[p_top_ray,t_top_ray,vpaflag_top_ray]=gx2cdf_ray([x_midtop x_top],w,k,lambda,m,s,'upper','mc_samples',2e3);
% time_ray=toc

% merge all parts
x=[x_bot x_midbot x_midtop x_top];
p_ncx2=[p_bot_ncx2 p_top_ncx2];
p_ray=[num2cell(p_bot_ray) p_top_ray];

% y-scaling:
% take log10 of all numeric and vpa numbers below the threshold,
% and scale above the threshold, for better visibility
p_ncx2_plot=log10(p_ncx2);
idx_below=p_ncx2<1e-3;
p_ncx2_plot(~idx_below)=1e3*p_ncx2(~idx_below)-4;

p_ray_plot=cellfun(@(x) double(log10(x)), p_ray);
idx_below=p_ray_plot<-3;
p_ray_plot(~idx_below)=1e3*[p_ray{~idx_below}]-4;

% x-scaling:

idx_trans=find(diff(idx_below));
x_bot=x(1:idx_trans(1)); % bottom part
x_mid=x(idx_trans(1)+1:idx_trans(2)); % middle part
x_top=x(idx_trans(2)+1:end); % top part

% take log10 of bottom:
x_bot_plot=log10(x_bot);
% 10x the middle and glue to the end of bottom:
x_mid_plot=x_bot_plot(end)+10*(x_mid-x_bot(end));
% 100* log10 of top, glue to end of mid:
x_top_plot=x_mid_plot(end)+100*(log10(x_top)-log10(x_mid(end)));
x_plot=[x_bot_plot x_mid_plot x_top_plot];

% separate double and vpa for plotting
first_vpa=find(~cellfun(@isnumeric, p_ray),1);

% plot
figure; hold on
plot(x_plot(first_vpa-1:end),p_ray_plot(first_vpa-1:end),'-ob','markersize',4,'markerfacecolor','w')
plot(x_plot(1:first_vpa-1),p_ray_plot(1:first_vpa-1),'-o','color',[0 .7 1],'markersize',4,'markerfacecolor','w')
plot(x_plot,p_ncx2_plot,'or','markersize',2,'markerfacecolor','r')
xline(x_mid_plot([1 length(x_midbot) end]))
yline(-3)
yline(log10(realmin))
xlim([-200 x_plot(end)])
ylim([-1.5e3 500])

set(gca,'xtick',[-200 -100 ...
    x_bot_plot(end)+10*([1 25]-x_bot(end))...
    x_mid_plot(end)+100*([2 3]-log10(x_mid(end)))...
    x_top_plot(end)])
set(gca,'xticklabel',{'10^{-200}','10^{-100}','1','25','10^2','10^3','10^4'})
set(gca,'ytick',[-1.5e3 -1e3 -500 -3 500])
set(gca,'yticklabel',{'10^{-1500}','10^{-1000}','10^{-500}','0.001','0.5'})

xlabel("$\chi'^2$",'interpreter','latex')
ylabel 'tail probability'
set(gca,'fontsize',13)

%% compare gx2 methods with ray method
% full gx2

w=[1 -5 2];
k=[1 2 3];
lambda=[2 3 7];
m=5;
s=10;

% find median and limits of middle part
med=gx2inv(0.5,w,k,lambda,m,s);
bot=gx2inv(0.01,w,k,lambda,m,s);
top=gx2inv(1-0.01,w,k,lambda,m,s);
[~,v]=gx2stat(w,k,lambda,m,s);

% bottom end
% x_tailbot=linspace(med-400*sqrt(v),bot,15);
x_bot=flip(-logspace(log10(abs(bot)),4.7,20));
x_midbot=linspace(bot,med,10);

% top end
x_midtop=linspace(med,top,10);
% x_top=linspace(top,med+200*sqrt(v),15);
x_top=logspace(log10(top),4.5,20);

% Davies method
tic
p_bot_gx2=gx2cdf([x_bot x_midbot],w,k,lambda,m,s,'method','imhof','AbsTol',0,'RelTol',1e-15);
p_top_gx2=gx2cdf([x_midtop x_top],w,k,lambda,m,s,'upper','method','imhof','AbsTol',0,'RelTol',1e-15);
time_gx2=toc

% ray method
tic
p_bot_ray=gx2cdf_ray([x_bot x_midbot],w,k,lambda,m,s,'mc_samples',2e3);
p_top_ray=gx2cdf_ray([x_midtop x_top],w,k,lambda,m,s,'upper','mc_samples',2e3);
time_ray=toc

% merge all parts
x=[x_bot x_midbot x_midtop x_top];
p_gx2=[p_bot_gx2 p_top_gx2];
p_ray=[p_bot_ray p_top_ray];

% y-scaling:
% take log10 of all numeric and vpa numbers below the threshold,
% and scale above the threshold, for better visibility
p_gx2_plot=log10(p_gx2);
idx_below=p_gx2<1e-2;
p_gx2_plot(~idx_below)=1e3*p_gx2(~idx_below)-12;

p_ray_plot=cellfun(@(x) double(log10(x)), p_ray);
p_ray_plot(~idx_below)=1e3*[p_ray{~idx_below}]-12;

% x-scaling:

idx_trans=find(diff(idx_below));
x_bot=x(1:idx_trans(1)); % bottom part
x_mid=x(idx_trans(1)+1:idx_trans(2)); % middle part
x_top=x(idx_trans(2)+1:end); % top part

% 70* log10 of bottom:
x_bot_plot=-70*log10(abs(x_bot));
% 1x the middle and glue to the end of bottom:
x_mid_plot=x_bot_plot(end)+1*(x_mid-x_bot(end));
% 70* log10 of top, glue to end of mid:
x_top_plot=x_mid_plot(end)+70*(log10(x_top)-log10(x_mid(end)));
x_plot=[x_bot_plot x_mid_plot x_top_plot];

% separate double and vpa
first_num=find(cellfun(@isnumeric, p_ray),1);
last_num=find(cellfun(@isnumeric, p_ray),1,'last');

% plot
figure; hold on
plot(x_plot(1:first_num),p_ray_plot(1:first_num),'-o','color',[0 .7 1],'markersize',4,'markerfacecolor','w')
plot(x_plot(last_num:end),p_ray_plot(last_num:end),'-o','color',[0 .7 1],'markersize',4,'markerfacecolor','w')
plot(x_plot(first_num:last_num),p_ray_plot(first_num:last_num),'-ob','markersize',4,'markerfacecolor','w')
plot(x_plot,p_gx2_plot,'or','markersize',2,'markerfacecolor','r')
med_plot=x_bot_plot(end)+1*(med-x_bot(end));
xline([x_mid_plot(1) med_plot x_mid_plot(end)])
yline([-2 log10(realmin)])
ylim([-1500 500])

xticks=[[-5 -4 -3]*70 ...
    x_bot_plot(end)+1*(x_mid([1 end])-x_bot(end))...
    x_mid_plot(end)+70*([3 4 5]-log10(x_mid(end)))];
set(gca,'xtick',xticks)
xlim(xticks([1 end]))
set(gca,'xticklabel',{'-10^{-5}','-10^{-4}','-10^{-3}','-64','56','10^3','10^4','10^5'})
set(gca,'ytick',[-1.5e3 -1e3 -500 -2 500])
set(gca,'yticklabel',{'10^{-1500}','10^{-1000}','10^{-500}','0.01','0.5'})

xlabel("$\tilde{\chi}$",'interpreter','latex')
ylabel 'tail probability'
set(gca,'fontsize',13)

%% illustration of Davies integrand
w=[1 -5 2];
k=[1 2 3];
lambda=[2 3 7];
m=5;
s=10;

x=-2e3;

% davies integrand
u=linspace(1e-3,.15,1e3);
f_davies=gx2cdf_imhof_integrand(u,x-m,w',k',lambda',s)/pi;

figure;
plot(u,f_davies,'-k')
axis([0 .15 -20 20])
set(gca,'xtick',[],'ytick',0,'fontsize',13)

%% 4d cubic function pdf using ray method

mu=[4;-2;3;2];
v=[  1 0 -1 0;
    0 8  4 0;
    -1 4  8 0;
    0 0  0 2];

fun_poly=@(x1,x2,x3,x4) x1.^3+x2.^2-x3.*x4;
fun_grad=@(x) [3*x(1)^2; 2*x(2); -x(4); -x(3)];

% norm_fun_inv(.01,mu,v,fun_poly)

x=linspace(-300,3000,50);
dx=x(2)-x(1);

tic
f_num=norm_fun_pdf(x,mu,v,fun_poly,'pdf_method','diff','fun_span',10,'fun_resol',100,'dx',1,'mc_samples',1e3);
time_num=toc

tic
f_ana=norm_fun_pdf(x,mu,v,fun_poly,'fun_grad',fun_grad,'fun_span',10,'fun_resol',100,'mc_samples',1e3);
time_ana=toc

% Monte-Carlo
tic
norm_mc=mvnrnd(mu,v,4e7);
fun_mc=arrayfun(fun_poly,norm_mc(:,1),norm_mc(:,2),norm_mc(:,3),norm_mc(:,4));
edges_mc=[x-dx/2 x(end)+dx/2];
f_mc=histcounts(fun_mc,edges_mc,'normalization','pdf');
time_mc=toc

figure; hold on
plot(x,f_ana,'-o','marker','.','markersize',4,'color',cobalt)
plot(x,f_num,'.','markersize',8,'color',sky)
plot(x,f_mc,'.','markersize',8,'color',orange)

set(gca,'yscale','log')
xlim([-300 2800])
ylim([1e-27 1e0])
xlabel("$f(\mbox{\boldmath $x$})$",'interpreter','latex')
ylabel 'pdf'
set(gca,'xtick',[0 1000 2000])
set(gca,'fontsize',13)


%% computing pdf: Imhof vs conv vs ray

w=[-2 1 1];
k=[2 1 2];
lambda=[7 5 3];
m=0;
s=3;

med=gx2inv(0.5,w,k,lambda,m,s);
[~,v]=gx2stat(w,k,lambda,m,s);

x_bot=-logspace(log10(4e3),log10(47),10);
x_mid=sort([linspace(-47,22,20) -7]);
x_top=logspace(log10(22),log10(2e3),10);
x=[x_bot x_mid x_top];

tic
f_conv=gx2pdf(x,w,k,lambda,m,s,'method','conv');
t_conv=toc

tic
f_im=gx2pdf_imhof(x,w,k,lambda,m,s,'AbsTol',0,'RelTol',1e-10);
t_im=toc

tic
f_ray=gx2pdf_ray(x,w,k,lambda,m,s,'mc_samples',1e5);
t_ray=toc

% plot
% scale x
x_bot_plot=20*sign(x_bot).*log10(abs(x_bot));
x_mid_plot=x_bot_plot(end)+x_mid-x_mid(1);
x_top_plot=x_mid_plot(end)+20*(log10(x_top)-log10(x_top(1)));
x_plot=[x_bot_plot x_mid_plot x_top_plot];

% scale y
mid_idx=f_ray>1e-3;

f_ray_plot=f_ray;
f_ray_plot(mid_idx)=f_ray_plot(mid_idx)*1e4-13;
f_ray_plot(~mid_idx)=log10(f_ray_plot(~mid_idx));

f_conv_plot=f_conv;
f_conv_plot(mid_idx)=f_conv_plot(mid_idx)*1e4-13;
f_conv_plot(~mid_idx)=log10(f_conv_plot(~mid_idx));

f_im_plot=f_im;
f_im_plot(mid_idx)=f_im_plot(mid_idx)*1e4-13;
f_im_plot(~mid_idx)=log10(f_im_plot(~mid_idx));

figure; hold on
plot(x_plot,f_ray_plot,'-o','markersize',6,'color',sky,'markerfacecolor','w','linewidth',1)
plot(x_plot,f_conv_plot,'ok','marker','.','markersize',5,'linewidth',1)
plot(x_plot,f_im_plot,'or','markersize',3,'color',orange,'linewidth',1)

yline([-3 log10(realmin)])
x_trans_plot=x_plot(find(diff(mid_idx)));
xline(x_trans_plot)
x_top_ticks=x_mid_plot(end)+20*(log10([1e3 1e4])-log10(x_top(1)));
xlim([-80 x_top_ticks(end)])
set(gca,'xtick',[-80 -60 x_trans_plot x_top_ticks])
set(gca,'XTickLabel',{'-10^4','-10^3','-47','22','10^3','10^4'})
ylim([-400 .04*1e4-13])
set(gca,'ytick',[-400 -200 -3 .04*1e4-13])
set(gca,'YTickLabel',{'10^{-400}','10^{-200}','0.001','0.04'})

xlabel("$\tilde{\chi}$",'interpreter','latex')
ylabel 'pdf'
set(gca,'fontsize',13)

%% extended d' check with vpa

mu_1=[0;0;0];

% both normals have the same covariance, so we can check against
% the true d' (Mahalanobis distance).
v=[1 .5 .7;
    .5  2  1 ;
    .7  1  3];

steps=logspace(0,5,20);
d_true=nan(size(steps));
d_gx2=nan(size(steps));
d_ray=nan(size(steps));

d_realmin=-2*norminv(realmin);

for i=1:length(steps)
    i
    mu_2=steps(i)*[1;1;1];

    results_gx2=classify_normals([mu_1,v],[mu_2,v],'method','gx2','plotmode',false);
    d_true(i)=results_gx2.norm_d_e;
    d_gx2(i)=results_gx2.norm_d_b;

    tic
    if d_true(i)<d_realmin
        results_ray=classify_normals([mu_1,v],[mu_2,v],'method','ray','AbsTol',0,'RelTol',1e-10,'plotmode',false);
    else
        results_ray=classify_normals([mu_1,v],[mu_2,v],'method','ray','vpa',true,'force_mc',true,'mc_samples',2e3,'plotmode',false);
    end
    toc
    d_ray(i)=results_ray.norm_d_b;
end

% for Monte Carlo
% steps_mc=linspace(1,13,10);
% d_true_mc=nan(size(steps_mc));
% d_mc=nan(size(steps_mc));
% 
% nsamp_mc=1e8;
% samp_1=mvnrnd(mu_1,v,nsamp_mc);
% samp_0=mvnrnd([0;0;0],v,nsamp_mc);
% 
% parfor i=1:length(steps_mc)
%     i
%     mu_2=steps_mc(i)*[1;1;1];
%     results_gx2=classify_normals([mu_1,v],[mu_2,v],'method','gx2','plotmode',0);
%     bd=results_gx2.norm_bd;
%     d_true_mc(i)=results_gx2.norm_d_e;
% 
%     samp_2=samp_0+mu_2';
%     d_mc(i)=-2*norminv(samp_value(samp_1,samp_2,bd,'vals',[0 1; 1 0])/(2*nsamp_mc));
% end

% calculate errors in estimates and plot
rel_err_gx2=abs(d_gx2-d_true)./d_true;
rel_err_gx2(~rel_err_gx2)=1e-18;
rel_err_ray=abs(d_ray-d_true)./d_true;
% rel_err_mc=abs(d_mc-d_true_mc)./d_true_mc;

figure; hold on
d_symmin=7e4;
% plot(d_true_mc,rel_err_mc,'-o','markersize',4,'color',.5*[1 1 1],'markerfacecolor',.5*[1 1 1])
% plot(d_true,rel_err_gx2,'-o','color',orange,'markersize',4,'markerfacecolor',orange)
first_vpa=find(d_true>d_realmin,1);
plot(d_true(first_vpa-1:end),rel_err_ray(first_vpa-1:end),'-o','color',cobalt,'markersize',4,'markerfacecolor',cobalt)
plot(d_true(1:first_vpa-1),rel_err_ray(1:first_vpa-1),'-o','color',sky,'markersize',4,'markerfacecolor',sky)
set(gca,'xscale','log')
set(gca,'yscale','log')
ylabel("rel. error $ | \hat{d}' - d' | / d' $",'Interpreter','latex');
xline([d_realmin d_symmin]) % largest computable d', corr. to the smallest possible error representable in double-precision
% yline(eps) % machine epsilon for double precision
xlabel("true $d'$",'interpreter','latex')
ylim([1e-18 1e-2])
set(gca,'ytick',[1e-18 1e-15 1e-10 1e-5])
set(gca,'yticklabel',{'0','10^{-15}','10^{-10}','10^{-5}'})
set(gca,'fontsize',13)