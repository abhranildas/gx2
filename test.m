
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

%% develop gx2pdf_ray

% ncx2
w=1;
k=4;
lambda=5;
m=0;
s=0;

x=linspace(0,50,1e3);
f_ray=gx2pdf_ray(x,w,k,lambda,m,s);

f_nc=ncx2pdf(x,k,lambda);

figure; hold on
plot(x,f_ray,'-k')
plot(x,f_nc)

% gx2

w=[1 -1];
k=[1 1];
lambda=[2 4];
m=0;
s=0;

x=linspace(-50,50,1e2);
tic
f_ray=gx2pdf_ray(x,w,k,lambda,m,s);
t_ray=toc;

tic
f_im=gx2pdf_imhof(x,w,k,lambda,m,s);
t_im=toc;

figure; hold on
plot(x,f_ray,'-k')
plot(x,f_im)


%% IFFT of characteristic function
% w=[-5 1 2];
% k=[1 2 1];
% lambda=[0 0 0];
% m=0;
% s=0;

w=1;
k=1;
lambda=0;
m=0;
s=0;

% w=0;
% k=1;
% lambda=0;
% m=0;
% s=1;

[f_conv,x_conv]=gx2pdf('full',w,k,lambda,m,s,'method','conv','xrange',150,'n_grid',1e6);
[f_ifft,x_ifft]=gx2pdf_ifft('full',w,k,lambda,m,s,'xrange',150,'n_grid',1e6);

figure;
hold on
plot(x_ifft,f_ifft)
plot(x_conv,f_conv)
set(gca,'yscale','log')
% xlim([-300 300])
% ylim([1e-20 1])

%% vpa for ray integration up to 3D

format long
dprime_true=100

mu_1=[0;0]; v=eye(2);
mu_2=dprime_true*[1;0];

results=classify_normals([mu_1,v],[mu_2,v],'AbsTol',0,'RelTol',1e-2,'vpa',true)

%% sym roots for SO
theta=-1; % theta = -1/0/1 give 2/1/0 valid roots respectively

q2=2*theta;
q1=theta+1;
q0=1;

r=valid_quad_roots(q2,q1,q0)

syms theta

q2=2*theta;
q1=theta+1;
q0=1;

r=valid_quad_roots(q2,q1,q0)

%% simplify psi(n)

format long
dprime_true=70

mu_1=[0;0;0]; v=eye(3);
mu_2=dprime_true*[1;0;0];

results=classify_normals([mu_1,v],[mu_2,v],'AbsTol',0,'RelTol',1e-2,'plotmode',0)
dprime=results.norm_d_b

