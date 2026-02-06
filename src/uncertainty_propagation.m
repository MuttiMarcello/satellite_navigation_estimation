% Spacecraft Guidance And Navigation
% Academic Year 2023/2024
% Assignment #2 exercise #1
% Mutti Marcello 220252 10698636

%% EX1
clear all; close all; clc;

% Set the default font size for axes labels, tick labels and legend
set(0,'DefaultAxesFontSize',20);
set(0,'DefaultLegendFontSize',20);

% Set the default linewidth
set(0,'DefaultLineLineWidth',2);

% Set the default text and legend interpreter
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% Necessary kernels loading
cspice_furnsh('assignment02.tm')

% Earth gravitational parameter
mu=cspice_bodvrd('Earth','GM',1);

% Data initialization
ets=cspice_str2et('2010-08-12 05:27:39.114 UTC');

% Mango mean state at ets
xs_M=[
    4622.232026629
    5399.3369588058
    -0.0212138165769957
    0.812221125483763
    -0.721512914578826
    7.42665302729053];

% Tango mean state at ets
xs_T=[
    4621.69343340281
    5399.26386352847
    -3.09039248714313
    0.813960847513811
    -0.719449862738607
    7.42706066911294];

% Covariance at ets
Ps=[
    +5.6e-7 +3.5e-7 -7.1e-8 0 0 0
    +3.5e-7 +9.7e-7 +7.6e-8 0 0 0
    -7.1e-8 +7.6e-8 +8.1e-8 0 0 0
    0 0 0 +2.8e-11 0 0
    0 0 0 0 +2.7e-11 0
    0 0 0 0 0 +9.6e-12];

% Initial STM conditions
Phi0=reshape(eye(6),[36,1]);

% Integrator settings
AbsTol=1e-12;   % Absolute tolerance
RelTol=1e-12;   % Relative tolerance
opt=odeset('AbsTol',AbsTol,'RelTol',RelTol);

% Mango semi-major axis and period
a1=mu/2*(mu/norm(xs_M(1:3))-norm(xs_M(4:6))^2/2)^(-1);
T1=2*pi*sqrt(a1^3/mu);
fprintf('Mango orbital period T1: %.4f s\n',T1)

% Tango semi-major axis and period
a2=mu/2*(mu/norm(xs_T(1:3))-norm(xs_T(4:6))^2/2)^(-1);
T2=2*pi*sqrt(a2^3/mu);
fprintf('Tango orbital period T2: %.4f s\n',T2)

% Number of orbits
N=10;

% Time grid
tt=linspace(ets,ets+N*T1,N+1);
tt_p=(0:N)*T1; % time grid for keplerian propagation

% LC
% state and STM propagation
[~,xx_M]=ode113(@(t,x) KEPSTM(t,x,mu),tt_p,[xs_M; Phi0],opt);
[~,xx_T]=ode113(@(t,x) KEPSTM(t,x,mu),tt_p,[xs_T; Phi0],opt);

% Means and covariances preallocation
xm_M_lc=zeros([length(tt),6]);
xm_T_lc=xm_M_lc;
PM_lc=zeros([6,6,length(tt)]);
PT_lc=PM_lc;

% 3sigmas allocations
sq9tr_M_lc=zeros(1,N+1);
sq9tr_M_lc(1)=3*sqrt(trace(Ps));
sq9tr_T_lc=sq9tr_M_lc;

for i=1:length(tt)
    xm_M_lc(i,:)=xx_M(i,1:6);   % Mango LC mean states
    xm_T_lc(i,:)=xx_T(i,1:6);   % Tango LC mean states

    % Mango LC covariances
    PM_lc(:,:,i)=reshape(xx_M(i,7:end),[6 6])*Ps*reshape(xx_M(i,7:end),[6 6])';

    % Tango LC covariances
    PT_lc(:,:,i)=reshape(xx_T(i,7:end),[6 6])*Ps*reshape(xx_T(i,7:end),[6 6])';
end

% UT
% covariances allocations
PM_ut=zeros([6,6,length(tt)]);
PM_ut(:,:,1)=Ps;
PT_ut=PM_ut;

% 3sigmas allocation
sq9tr_M_ut=sq9tr_M_lc;
sq9tr_T_ut=sq9tr_M_ut;

% UT tuning parameters
n=length(xs_M);
alpha=1e-1;
beta=2;
lambda=alpha^2*n-n;

% UT weights
Wm=zeros(2*n+1,1);
Wc=Wm;
Wm(1)=lambda/(lambda+n);
Wm(2:end)=1/(2*(lambda+n));
Wc(1)=Wm(1)+(1-alpha^2+beta);
Wc(2:end)=Wm(2:end);

XC=sqrtm((n+lambda)*Ps);

% sigma point matrix allocation
YM=zeros(N+1,6,2*n+1);
YT=YM;
for i=1:2*n+1

    % sigma point computation
    if i==1
        XM=xs_M;
        XT=xs_T;
    elseif i<=n+1
        XM=xs_M+XC(:,i-1);
        XT=xs_T+XC(:,i-1);
    else
        XM=xs_M-XC(:,i-1-n);
        XT=xs_T-XC(:,i-1-n);
    end

    % sigma point propagation and collection
    [~,Y]=ode113(@(t,x) KEP(t,x,mu),tt_p,XM,opt);
    YM(:,:,i)=Y;
    [~,Y]=ode113(@(t,x) KEP(t,x,mu),tt_p,XT,opt);
    YT(:,:,i)=Y;
end

% Tango UT mean states
xm_M_ut = sum(YM.*reshape(Wm,1,1,[]),3);
xm_M_ut(1,:)=xs_M;

% Mango UT mean states
xm_T_ut = sum(YT.*reshape(Wm,1,1,[]),3);
xm_T_ut(1,:)=xs_T;


for i=2:N+1
    for j=1:2*n+1

        % Mango UT covariances
        PM_ut(:,:,i)=PM_ut(:,:,i)+Wc(j)*(xm_M_ut(i,:)-YM(i,:,j))'*(xm_M_ut(i,:)-YM(i,:,j));
        
        % Tango UT covariance
        PT_ut(:,:,i)=PT_ut(:,:,i)+Wc(j)*(xm_T_ut(i,:)-YT(i,:,j))'*(xm_T_ut(i,:)-YT(i,:,j));
    end

    % LC and UT 3sigma computation
    sq9tr_M_lc(i)=3*sqrt(trace(PM_lc(:,:,i)));
    sq9tr_M_ut(i)=3*sqrt(trace(PM_ut(:,:,i)));
    sq9tr_T_lc(i)=3*sqrt(trace(PT_lc(:,:,i)));
    sq9tr_T_ut(i)=3*sqrt(trace(PT_ut(:,:,i)));
end


% Methods propagated position and velocity comparisons
figure
subplot(2,2,1)
plot(0:N,sqrt(sum((xm_M_lc(:,1:3)-xm_M_ut(:,1:3)).^2,2)),'o-')
grid minor
xlim([0 N])
ylabel('$\left|\hat{\mathbf{r}}^{LC}-\hat{\mathbf{r}}^{UT}\right|\,[km]$')
title('Mango')
subplot(2,2,3)
plot(0:N,sqrt(sum((xm_M_lc(:,4:6)-xm_M_ut(:,4:6)).^2,2)),'o-')
grid minor
xlim([0 N])
ylabel('$\left|\hat{\mathbf{v}}^{LC}-\hat{\mathbf{v}}^{UT}\right|\,[km/s]$')
xlabel('$n\,[-]$')
subplot(2,2,2)
plot(0:N,sqrt(sum((xm_T_lc(:,1:3)-xm_T_ut(:,1:3)).^2,2)),'o-')
grid minor
xlim([0 N])
title('Tango')
subplot(2,2,4)
plot(0:N,sqrt(sum((xm_T_lc(:,4:6)-xm_T_ut(:,4:6)).^2,2)),'o-')
grid minor
xlim([0 N])
xlabel('$n\,[-]$')

% Methods propagated covariances comparisons
figure
subplot(2,2,1)
plot(0:N,sq9tr_M_lc,'bo-')
hold on
plot(0:N,sq9tr_M_ut,'ro-')
xlim([0 N])
grid minor
ylabel('$3\sqrt{tr\left(P\right)}$')
legend('$3\sqrt{tr\left(P_M^{LC}\right)}$','$3\sqrt{tr\left(P_M^{UT}\right)}$','fontsize',23,'location','northwest')
title('Mango')
subplot(2,2,2)
plot(0:N,sq9tr_T_lc,'bo-')
hold on
plot(0:N,sq9tr_T_ut,'ro-')
xlim([0 N])
grid minor
legend('$3\sqrt{tr\left(P_T^{LC}\right)}$','$3\sqrt{tr\left(P_T^{UT}\right)}$','fontsize',23,'location','northwest')
title('Tango')
subplot(2,2,3)
plot(0:N,abs(sq9tr_M_lc-sq9tr_M_ut),'o-')
xlim([0 N])
grid minor
ylabel('$\left|3\sqrt{tr\left(P^{LC}\right)}-3\sqrt{tr\left(P^{UT}\right)}\right|$')
xlabel('$n\,[-]$')
subplot(2,2,4)
plot(0:N,abs(sq9tr_T_lc-sq9tr_T_ut),'o-')
xlim([0 N])
grid minor
xlabel('$n\,[-]$')


%% EX2
close all; clc;

% LC computed relative position
r_rel_lc=sqrt(sum((xm_M_lc(:,1:3)-xm_T_lc(:,1:3)).^2,2));

% UT computed relative position
r_rel_ut=sqrt(sum((xm_M_ut(:,1:3)-xm_T_ut(:,1:3)).^2,2));

% LC computed Psum
P_sum_lc=PM_lc+PT_lc;

% UT computed Psum
P_sum_ut=PM_ut+PT_ut;

% Eigenvalues allocation
eig_lc=zeros(length(tt),1);
eig_ut=eig_lc;

for i=1:N+1
    
    % LC covariances max eig extraction
    eig_lc(i)=max(eig(P_sum_lc(:,:,i)));

    % UT covariances max eig extraction
    eig_ut(i)=max(eig(P_sum_ut(:,:,i)));
end

% LC critical condition
N_cr_lc=find(r_rel_lc<3*sqrt(eig_lc),1)-1;

% UT critical condition
N_cr_ut=find(r_rel_ut<3*sqrt(eig_ut),1)-1;

% Methods critical warning comparison
figure
subplot(1,2,1)
plot(0:N,r_rel_lc,'-o',0:N,3*sqrt(eig_lc),'-o')
hold on
xline(N_cr_lc,'k-.','LineWidth',0.5)
grid minor
xlabel('$n\,[-]$')
ylabel('$\Delta r ,\quad 3\sqrt{max\left(\lambda_i(P_{sum})\right)}\quad [km]$')
legend('$\Delta r$','$3\sqrt{max\left(\lambda_i(P_{sum})\right)}$','$n_c$','fontsize',25,'location','northeast')
title('LC')
subplot(1,2,2)
plot(0:N,r_rel_ut,'-o',0:N,3*sqrt(eig_ut),'-o')
hold on
xline(N_cr_ut,'k-.','LineWidth',0.5)
grid minor
xlabel('$n\,[-]$')
ylabel('$\Delta r ,\quad 3\sqrt{max\left(\lambda_i(P_{sum})\right)}\quad [km]$')
legend('$\Delta r$','$3\sqrt{max\left(\lambda_i(P_{sum})\right)}$','$n_c$','fontsize',25,'location','northeast')
title('UT')

% Methods Psum comparison
figure
plot(0:N,(3*sqrt(eig_ut)-3*sqrt(eig_lc)),'-o')
grid minor
xlabel('$n\,[-]$')
ylabel('$3\sqrt{max\left(\lambda_i(P_{sum}^{UT})\right)}-3\sqrt{max\left(\lambda_i(P_{sum}^{LC})\right)}$')
title('LC vs UT')

%% EX3
close all; clc;

% Number of monte carlo samples
n_mc=250;

% MC mean states allocation
XM_mc=zeros(N+1,6,n_mc);
XT_mc=XM_mc;

% MC covariances allocation
PM_mc=zeros(6,6,N+1);
PT_mc=PM_mc;    

% Initial random samples generation
for i=1:n_mc
    % Mango initial random samples
    x0_M_mc=mvnrnd(xs_M,Ps)';
    [~,xx]=ode113(@(t,x) KEP(t,x,mu),tt_p,x0_M_mc,opt);
    XM_mc(:,:,i)=xx;

    % Tango initial random samples
    x0_T_mc=mvnrnd(xs_T,Ps)';
    [~,xx]=ode113(@(t,x) KEP(t,x,mu),tt_p,x0_T_mc,opt);
    XT_mc(:,:,i)=xx;
end

% MC sample mean states
xm_M_mc=mean(XM_mc,3);
xm_T_mc=mean(XT_mc,3);

% MC sample covariances
for i=1:N+1
    for j=1:n_mc
        PM_mc(:,:,i)=PM_mc(:,:,i)+(xm_M_mc(i,:)-XM_mc(i,:,j))'*(xm_M_mc(i,:)-XM_mc(i,:,j));
        PT_mc(:,:,i)=PT_mc(:,:,i)+(xm_T_mc(i,:)-XT_mc(i,:,j))'*(xm_T_mc(i,:)-XT_mc(i,:,j));
    end
end

% Mango MC sample covariance
PM_mc=PM_mc./(n_mc-1);

% Mango MC sample covariance
PT_mc=PT_mc./(n_mc-1);
 
% position and velocity max eig method allocation
lr_lc=zeros(2,N+1); % position
lr_ut=lr_lc;
lr_mc=lr_lc;
lv_lc=lr_lc;        % velocity
lv_ut=lr_lc;
lv_mc=lr_lc;

% max eig computation per method
for i=1:N+1
    lr_lc(1,i)=3*sqrt(max(eig(PM_lc(1:3,1:3,i))));
    lr_lc(2,i)=3*sqrt(max(eig(PT_lc(1:3,1:3,i))));

    lr_ut(1,i)=3*sqrt(max(eig(PM_ut(1:3,1:3,i))));
    lr_ut(2,i)=3*sqrt(max(eig(PT_ut(1:3,1:3,i))));

    lr_mc(1,i)=3*sqrt(max(eig(PM_mc(1:3,1:3,i))));
    lr_mc(2,i)=3*sqrt(max(eig(PT_mc(1:3,1:3,i))));

    lv_lc(1,i)=3*sqrt(max(eig(PM_lc(4:6,4:6,i))));
    lv_lc(2,i)=3*sqrt(max(eig(PT_lc(4:6,4:6,i))));

    lv_ut(1,i)=3*sqrt(max(eig(PM_ut(4:6,4:6,i))));
    lv_ut(2,i)=3*sqrt(max(eig(PT_ut(4:6,4:6,i))));

    lv_mc(1,i)=3*sqrt(max(eig(PM_mc(4:6,4:6,i))));
    lv_mc(2,i)=3*sqrt(max(eig(PT_mc(4:6,4:6,i))));
end

% Plotted revolutions identifiers
id=[1 4 7 10]+1;

% LVLH frame rotation matrix allocations
RM=zeros(3,3);
RT=RM;

% LVLH Mango methods covariances allocation
PM_lc_lvlh=zeros(3,3,length(id));
PM_mc_lvlh=PM_lc_lvlh;
PM_ut_lvlh=PM_lc_lvlh;


% LVLH Tango methods covariances allocation
PT_lc_lvlh=PM_lc_lvlh;
PT_mc_lvlh=PM_lc_lvlh;
PT_ut_lvlh=PM_lc_lvlh;

% LVLH Mango mean states allocation
rm_M_lc_lvlh=zeros(3,length(id));
rm_M_mc_lvlh=zeros(3,length(id));
rm_M_ut_lvlh=zeros(3,length(id));

% LVLH Tango mean states allocation
rm_T_lc_lvlh=zeros(3,length(id));
rm_T_mc_lvlh=zeros(3,length(id));
rm_T_ut_lvlh=zeros(3,length(id));

% LVLH Mango MC samples allocation
XM_mc_lvlh=zeros(3,n_mc,length(id));

% LVLH Tango MC samples allocation
XT_mc_lvlh=XM_mc_lvlh;

for i=1:length(id)

    % Mango LVLH rotation matrix
    RM(1,:)=xm_M_lc(id(i),1:3)/norm(xm_M_lc(id(i),1:3));
    RM(3,:)=cross(RM(1,:),xm_M_lc(id(i),4:6)/norm(xm_M_lc(id(i),4:6)))/norm(cross(RM(1,:),xm_M_lc(id(i),4:6)/norm(xm_M_lc(id(i),4:6))));
    RM(2,:)=cross(RM(3,:),RM(1,:));

    % Tango LVLH rotation matrix
    RT(1,:)=xm_T_lc(id(i),1:3)/norm(xm_T_lc(id(i),1:3));
    RT(3,:)=cross(RT(1,:),xm_T_lc(id(i),4:6)/norm(xm_T_lc(id(i),4:6)))/norm(cross(RT(1,:),xm_T_lc(id(i),4:6)/norm(xm_T_lc(id(i),4:6))));
    RT(2,:)=cross(RT(3,:),RT(1,:));

    % Mango LVLH mean states computation
    rm_M_lc_lvlh(:,i)=RM*(xm_M_lc(id(i),1:3)-xm_M_lc(id(i),1:3))';
    rm_M_mc_lvlh(:,i)=RM*(xm_M_mc(id(i),1:3)-xm_M_lc(id(i),1:3))';
    rm_M_ut_lvlh(:,i)=RM*(xm_M_ut(id(i),1:3)-xm_M_lc(id(i),1:3))';
    
    % Tango LVLH mean states computation
    rm_T_lc_lvlh(:,i)=RT*(xm_T_lc(id(i),1:3)-xm_T_lc(id(i),1:3))';
    rm_T_mc_lvlh(:,i)=RT*(xm_T_mc(id(i),1:3)-xm_T_lc(id(i),1:3))';
    rm_T_ut_lvlh(:,i)=RT*(xm_T_ut(id(i),1:3)-xm_T_lc(id(i),1:3))';

    % Mango LVLH covariances computation
    PM_lc_lvlh(:,:,i)=RM*reshape(PM_lc(1:3,1:3,id(i)),[3 3])*RM';
    PM_mc_lvlh(:,:,i)=RM*reshape(PM_mc(1:3,1:3,id(i)),[3 3])*RM';
    PM_ut_lvlh(:,:,i)=RM*reshape(PM_ut(1:3,1:3,id(i)),[3 3])*RM';

    % Tango LVLH covariances computation
    PT_lc_lvlh(:,:,i)=RT*reshape(PT_lc(1:3,1:3,id(i)),[3 3])*RT';
    PT_mc_lvlh(:,:,i)=RT*reshape(PT_mc(1:3,1:3,id(i)),[3 3])*RT';
    PT_ut_lvlh(:,:,i)=RT*reshape(PT_ut(1:3,1:3,id(i)),[3 3])*RT';

    % Mango LVLH MC samples computation
    XM_mc_lvlh(:,:,i)=RM*(reshape(XM_mc(id(i),1:3,:),[3 n_mc])-xm_M_lc(id(i),1:3)');
    
    % Tango LVLH MC samples computation
    XT_mc_lvlh(:,:,i)=RT*(reshape(XT_mc(id(i),1:3,:),[3 n_mc])-xm_T_lc(id(i),1:3)');
end

% Selected revolutions covariance ellipses rapresentation
for i=1:length(id)

    % Mango
    figure(10*i)

    % MC samples
    for j=1:n_mc    
        plot(XM_mc_lvlh(1,j,i),XM_mc_lvlh(2,j,i),'x','Color',[.7 .7 .7],'LineWidth',0.5);
        hold on
    end
    plot(rm_M_mc_lvlh(1,i),rm_M_mc_lvlh(2,i),'bo')
    grid minor
    plot(rm_M_lc_lvlh(1,i),rm_M_lc_lvlh(2,i),'rd')
    plot(rm_M_ut_lvlh(1,i),rm_M_ut_lvlh(2,i),'gs')

    % MC 3sigma ellipse
    xy=cov_ellipse(PM_mc_lvlh(1:3,1:3,i));
    mcM=plot(xy(1,:)+rm_M_mc_lvlh(1,i),xy(2,:)+rm_M_mc_lvlh(2,i),'b','DisplayName','MC');

    % LC 3sigma ellipse
    xy=cov_ellipse(PM_lc_lvlh(1:3,1:3,i));
    lcM=plot(xy(1,:)+rm_M_lc_lvlh(1,i),xy(2,:)+rm_M_lc_lvlh(2,i),'r','DisplayName','LC');

    % UT 3sigma ellipse
    xy=cov_ellipse(PM_ut_lvlh(1:3,1:3,i));
    utM=plot(xy(1,:)+rm_M_ut_lvlh(1,i),xy(2,:)+rm_M_ut_lvlh(2,i),'g','DisplayName','UT');

    xline(0,'-.','LineWidth',0.5)
    yline(0,'-.','LineWidth',0.5)
    xlim([-0.01 0.01])
    ylim([-2 2])
    xlabel('$\hat{\mathbf{i}}_{LVLH}\quad [km]$')
    ylabel('$\hat{\mathbf{j}}_{LVLH}\quad [km]$')
    tit=sprintf('Mango, N=%.0f',id(i)-1);
    legend([mcM lcM utM])
    title(tit)

    % Tango
    figure(100*i)

    % MC samples
    for j=1:n_mc    
        plot(XT_mc_lvlh(1,j,i),XT_mc_lvlh(2,j,i),'x','Color',[.7 .7 .7],'LineWidth',0.5)
        hold on
    end
    plot(rm_T_mc_lvlh(1,i),rm_T_mc_lvlh(2,i),'bo')
    grid minor
    plot(rm_T_lc_lvlh(1,i),rm_T_lc_lvlh(2,i),'rd')
    plot(rm_T_ut_lvlh(1,i),rm_T_ut_lvlh(2,i),'gs')

    % MC 3sigma ellipses
    xy=cov_ellipse(PT_mc_lvlh(1:3,1:3,i));
    mcT=plot(xy(1,:)+rm_T_mc_lvlh(1,i),xy(2,:)+rm_T_mc_lvlh(2,i),'b','DisplayName','MC');

    % LC 3sigma ellipses
    xy=cov_ellipse(PT_lc_lvlh(1:3,1:3,i));
    lcT=plot(xy(1,:)+rm_T_lc_lvlh(1,i),xy(2,:)+rm_T_lc_lvlh(2,i),'r','DisplayName','LC');

    % UT 3sigma ellipses
    xy=cov_ellipse(PT_ut_lvlh(1:3,1:3,i));
    utT=plot(xy(1,:)+rm_T_ut_lvlh(1,i),xy(2,:)+rm_T_ut_lvlh(2,i),'g','DisplayName','UT');
    
    xline(0,'-.','LineWidth',0.5)
    yline(0,'-.','LineWidth',0.5)
    xlim([-0.01 0.01])
    ylim([-2 2])
    xlabel('$\hat{\mathbf{i}}_{LVLH}\quad [km]$')
    ylabel('$\hat{\mathbf{j}}_{LVLH}\quad [km]$')
    tit=sprintf('Tango, N=%.0f',id(i)-1);
    legend([mcT lcT utT])
    title(tit)
end

% Mango & Tango position covariance evolution
figure
subplot(1,2,1)
plot(0:N,lr_lc(1,:),'o-',0:N,lr_ut(1,:),'o--',0:N,lr_mc(1,:),'o:')
grid minor
xlabel('$n\,[-]$')
ylabel('$3\sqrt{max\left(\lambda_i(P_{r,M})\right)}\quad [km]$')
legend('LC','UT','MC','location','southeast')
title('Mango')
axes('Position',[0.19,0.7,0.1,0.2])
box on;
plot(N-1:N,lr_lc(1,end-1:end),'o-',N-1:N,lr_ut(1,end-1:end),'o--',N-1:N,lr_mc(1,end-1:end),'o:')
ylim([min([lr_lc(1,end),lr_ut(1,end),lr_mc(1,end)]) max([lr_lc(1,end),lr_ut(1,end),lr_mc(1,end)])])
grid minor
fprintf('Mango, maximum 3sigma position displacement wrt MC: %.2f %%\n',100*abs(max(lr_mc(1,end)-lr_ut(1,end),lr_mc(1,end)-lr_lc(1,end)))/lr_mc(1,end))
subplot(1,2,2)
plot(0:N,lr_lc(2,:),'o-',0:N,lr_ut(2,:),'o--',0:N,lr_mc(2,:),'o:')
grid minor
xlabel('$n\,[-]$')
ylabel('$3\sqrt{max\left(\lambda_i(P_{r,T})\right)}\quad [km]$')
legend('LC','UT','MC','location','southeast')
title('Tango')
axes('Position',[0.63,0.7,0.1,0.2])
box on;
plot(N-1:N,lr_lc(2,end-1:end),'o-',N-1:N,lr_ut(2,end-1:end),'o--',N-1:N,lr_mc(2,end-1:end),'o:')
ylim([min([lr_lc(2,end),lr_ut(2,end),lr_mc(2,end)]) max([lr_lc(2,end),lr_ut(2,end),lr_mc(2,end)])])
grid minor
fprintf('Tango, maximum 3sigma position displacement wrt MC: %.2f %%\n',100*abs(max(lr_mc(2,end)-lr_ut(2,end),lr_mc(2,end)-lr_lc(2,end)))/lr_mc(2,end))


% Mango & Tango velocity covariance evolution
figure
subplot(1,2,1)
plot(0:N,lv_lc(1,:),'o-',0:N,lv_ut(1,:),'o--',0:N,lv_mc(1,:),'o:')
grid minor
xlabel('$n\,[-]$')
ylabel('$3\sqrt{max\left(\lambda_i(P_{v,M})\right)}\quad [km/s]$')
legend('LC','UT','MC','location','southeast')
title('Mango')
axes('Position',[0.19,0.67,0.1,0.2])
box on;
plot(N-1:N,lv_lc(1,end-1:end),'o-',N-1:N,lv_ut(1,end-1:end),'o--',N-1:N,lv_mc(1,end-1:end),'o:')
ylim([min([lv_lc(1,end),lv_ut(1,end),lv_mc(1,end)]) max([lv_lc(1,end),lv_ut(1,end),lv_mc(1,end)])])
grid minor
fprintf('Mango, maximum 3sigma velocity displacement wrt MC: %.2f %%\n',100*abs(max(lv_mc(1,end)-lv_ut(1,end),lv_mc(1,end)-lv_lc(1,end)))/lv_mc(1,end))
subplot(1,2,2)
plot(0:N,lv_lc(2,:),'o-',0:N,lv_ut(2,:),'o--',0:N,lv_mc(2,:),'o:')
grid minor
xlabel('$n\,[-]$')
ylabel('$3\sqrt{max\left(\lambda_i(P_{v,T})\right)}\quad [km/s]$')
legend('LC','UT','MC','location','southeast')
title('Tango')
axes('Position',[0.63,0.67,0.1,0.2])
box on;
plot(N-1:N,lv_lc(2,end-1:end),'o-',N-1:N,lv_ut(2,end-1:end),'o--',N-1:N,lv_mc(2,end-1:end),'o:')
ylim([min([lv_lc(2,end),lv_ut(2,end),lv_mc(2,end)]) max([lv_lc(2,end),lv_ut(2,end),lv_mc(2,end)])])
grid minor
fprintf('Tango, maximum 3sigma velocity displacement wrt MC: %.2f %%\n',100*abs(max(lv_mc(2,end)-lv_ut(2,end),lv_mc(2,end)-lv_lc(2,end)))/lv_mc(2,end))

%% User plot setting removal
set(0,'DefaultAxesFontSize','remove');
set(0,'DefaultLegendFontSize','remove');
set(0,'DefaultLineLineWidth','remove');
set(groot,'defaultTextInterpreter','remove');
set(groot,'defaultLegendInterpreter','remove');
cspice_kclear;

%% Functions

function [rhs] = KEPSTM(~,x,mu)
%     Computes keplerian equations of motion RHS together with relative STM
%     propagation RHS
%     Example: [rhs] = KEPSTM(t, x, mu)
%     INPUTS:
%         t  [1x1]  reference epoch (omissible)
%         x  [1x42] state and state transition matrix
%         mu [1x1]  KEPSTM gravitational constant
%     OUTPUTS:
%         rhs [42x1] equations of motion and STM propagation RHS

    % State and STM verticality check
    [m,n]=size(x);
    if n>m
        x=x';
    end

    % State Variables Extraction
    rx=x(1);
    ry=x(2);
    rz=x(3);
    vx=x(4);
    vy=x(5);
    vz=x(6);

    rr=[rx ry rz].';
    
    % STM Extraction and Reshaping 
    Phi=reshape(x(7:end),6,6);
    
    Ar=[[(3*mu*rx*abs(rx)*(rx + (rx)))/(2*(rx*(rx))^(1/2)*(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(5/2)) - mu/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(3/2),                                                (3*mu*rx*abs(ry)*(ry + (ry)))/(2*(ry*(ry))^(1/2)*(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(5/2)),                                                (3*mu*rx*abs(rz)*(rz + (rz)))/(2*(rz*(rz))^(1/2)*(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(5/2))]
       [                                               (3*mu*ry*abs(rx)*(rx + (rx)))/(2*(rx*(rx))^(1/2)*(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(5/2)), (3*mu*ry*abs(ry)*(ry + (ry)))/(2*(ry*(ry))^(1/2)*(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(5/2)) - mu/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(3/2),                                                (3*mu*ry*abs(rz)*(rz + (rz)))/(2*(rz*(rz))^(1/2)*(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(5/2))]
       [                                               (3*mu*rz*abs(rx)*(rx + (rx)))/(2*(rx*(rx))^(1/2)*(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(5/2)),                                                (3*mu*rz*abs(ry)*(ry + (ry)))/(2*(ry*(ry))^(1/2)*(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(5/2)), (3*mu*rz*abs(rz)*(rz + (rz)))/(2*(rz*(rz))^(1/2)*(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(5/2)) - mu/(abs(rx)^2 + abs(ry)^2 + abs(rz)^2)^(3/2)]];

    A=zeros(6,6);
    A(1:3,4:6)=eye(3);
    A(4:6,1:3)=Ar;

    % STM rhs computation
    Phidot=A*Phi;

    % RHS assembly
    rhs=zeros(42,1);
    rhs(1:3)=[vx vy vz];
    rhs(4:6)=-mu.*rr./norm(rr)^3;
    rhs(7:end)=reshape(Phidot,[36,1]);
end

function [rhs] = KEP(~,x,mu)
%     Computes keplerian equations of motion RHS propagation RHS
%     Example: [rhs] = KEP(t, x, mu)
%     INPUTS:
%         t  [1x1]  reference epoch (omissible)
%         x  [1x6]  state
%         mu [1x1]  KEPSTM gravitational constant
%     OUTPUTS:
%         rhs [6x1] equations of motion and STM propagation RHS

    % State verticality check
    [m,n]=size(x);
    if n>m
        x=x';
    end

    % State Variables Extraction
    rx=x(1);
    ry=x(2);
    rz=x(3);
    vx=x(4);
    vy=x(5);
    vz=x(6);

    rr=[rx ry rz].';
    
    % RHS assembly
    rhs=zeros(6,1);
    rhs(1:3)=[vx vy vz];
    rhs(4:6)=-mu.*rr./norm(rr)^3;
end

function [xye] = cov_ellipse(P) 
%     Computes 3 sigma covariance ellipse of 2x2 northwest minor of P
%     Example: [xye] = cov_ellipse(P)
%     INPUTS:
%         P   [nxn]  covariance n>1
%     OUTPUTS:
%         xye [2x50] 3sigma covariance ellipse contour coordinates
    
    xye=zeros(2,50);
    t=linspace(0,2*pi,50);

    % Eigenvalues and vectors extraction
    [V,L]=eig(P(1:2,1:2));

    % square of 3sigma axes
    L1=9*abs(L(1,1));
    L2=9*abs(L(2,2));
    V1=V(:,1);

    % ellipse tilt angle
    th=atan2(V1(2),V1(1));

    % coordinates
    xye(1,:)=sqrt(L1)*cos(th)*cos(t)-sqrt(L2)*sin(th)*sin(t);
    xye(2,:)=sqrt(L1)*sin(th)*cos(t)+sqrt(L2)*cos(th)*sin(t);
end
