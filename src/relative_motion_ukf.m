% Spacecraft Guidance And Navigation
% Academic Year 2023/2024
% Assignment #2 exercise #3
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

% Necessary kernels and SPG4 loading
cspice_furnsh('assignment02.tm')
addpath('sgp4');

% Earth gravitational parameter
mu=cspice_bodvrd('Earth','GM',1);

% Data inizialization
ets=cspice_str2et('2010-08-12 05:27:39.114 UTC');
et0=cspice_str2et('2010-08-12 05:30:00.000 UTC');
etf=cspice_str2et('2010-08-12 06:30:00.000 UTC');

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

% Integrator settings
AbsTol=1e-12;   % Absolute tolerance
RelTol=1e-12;   % Relative tolerance
opt=odeset('AbsTol',AbsTol,'RelTol',RelTol);

% Mango state at et0
[~,xx0]=ode113(@(t,x) KEP(t,x,mu),[ets et0],xs_M,opt);
x0_M=xx0(end,:)';

% Time grid
tm=5;   % sampling time
tt=et0:tm:etf;

% keplerian Mango propagated states
[~,xx_M]=ode113(@(t,x) KEP(t,x,mu),tt,x0_M,opt);

% svalbard
str_S='SVALBARD';
min_el_S=5; % [deg] minimum elevation requirement
R_S=diag([deg2rad(125*1e-3) deg2rad(125*1e-3) 0.01].^2); % radar noise matrix
% anticipated keplerian measurements
[az_M_S, el_M_S, r_M_S]=antenna_pointing(str_S,tt,xx_M);

% visibility identifiers
S_vis_id_M=el_M_S>deg2rad(min_el_S);

% anticipated visible keplerian measurements
meas_kep_S_M=[az_M_S(S_vis_id_M==1);el_M_S(S_vis_id_M==1);r_M_S(S_vis_id_M==1)];

% anticipated visibility epochs
tt_S_M=tt(S_vis_id_M==1);

% allocated svalbard visibility window
VW_S_M=[];
VW_S_M(1,:)=tt(strfind(S_vis_id_M,[0 1])+1);  % first row opening of window
VW_S_M(2,:)=tt(strfind(S_vis_id_M,[1 0]));    % second row closing of window

fprintf('Mango\n')
for i=1:size(VW_S_M,2)
    fprintf('%s expected visibility window %.0f from %s to %s\n',str_S,i,cspice_et2utc(VW_S_M(1,i),'C',3),cspice_et2utc(VW_S_M(2,i),'C',3))
end
fprintf('Total %s expected visibility: %.2f min, %.0f measurements\n\n',str_S,sum(VW_S_M(2,:)-VW_S_M(1,:))/60,sum(S_vis_id_M))

typerun='u';  % user-provided inputs to SGP4 Matlab function
opsmode='a';  % afspc approach ('air force space command')
whichconst=72;  % WGS72 constants (radius, gravitational parameter)

% Mango TLE and satrec
M_longstr1='1 36599U 10028B   10224.22752732 -.00000576  00000-0 -16475-3 0  9998';
M_longstr2='2 36599 098.2803 049.5758 0043871 021.7908 338.5082 14.40871350  8293';
M_rec=twoline2rv(M_longstr1, M_longstr2, typerun,'e', opsmode, whichconst);

[year,mon,day,hr,min,sec]=invjday(M_rec.jdsatepoch,M_rec.jdsatepochf);
et_M_TLE_str=sprintf('%d-%02d-%02d %02d:%02d:%02.6f',[year,mon,day,hr,min,sec]);

% Mango reference TLE epoch
et_M_TLE=cspice_str2et(et_M_TLE_str);

fprintf('Mango satellite ID: %d\n',M_rec.satnum);
fprintf('TLE reference epoch: UTC %s\n\n',et_M_TLE_str);

% svalbard measurements
[xx_ref_M,meas_S_M,tt_meas_S_M,S_vis_id_M]=station_measurements(M_rec,str_S,tt_S_M,et_M_TLE);

% station apparent motion
figure
plot(rad2deg(meas_kep_S_M(1,:)),rad2deg(meas_kep_S_M(2,:)),'o')
hold on
plot(rad2deg(meas_S_M(1,:)),rad2deg(meas_S_M(2,:)),'x')
yline(min_el_S,'-.','LineWidth',0.5)
grid minor
xlim([10 145])
ylim([0 round(max(rad2deg([meas_kep_S_M(2,:), meas_S_M(2,:)])),-1)+5*(round(max(rad2deg([meas_kep_S_M(2,:), meas_S_M(2,:)])),-1)<max(rad2deg([meas_kep_S_M(2,:), meas_S_M(2,:)])))])
xlabel('$Az\, [deg]$')
ylabel('$El\, [deg]$')
legend('Expected Keplerian Measurements','Measurements','Location','south','fontsize',25)
title(sprintf('Mango - %s',str_S))

% station measurements
figure
subplot(3,1,1)
plot(5*(0:(length(tt_meas_S_M)-1))/60,rad2deg(meas_S_M(1,:)),'rx')
grid minor
xlim([0 5*(length(tt_meas_S_M)-1)/60])
ylabel('$Az\, [deg]$')
title(str_S)
subplot(3,1,2)
plot(5*(0:(length(tt_meas_S_M)-1))/60,rad2deg(meas_S_M(2,:)),'rx')
grid minor
xlim([0 5*(length(tt_meas_S_M)-1)/60])
ylabel('$El\, [deg]$')
subplot(3,1,3)
plot(5*(0:(length(tt_meas_S_M)-1))/60,(meas_S_M(3,:)),'rx')
grid minor
xlim([0 5*(length(tt_meas_S_M)-1)/60])
ylabel('$r\, [km]$')
% xlabel('Measurements\, [-]')
xlabel(sprintf('Minutes past %s',cspice_et2utc(tt_meas_S_M(1),'C',0)))

% UT tuning parameters
n=length(xs_M);
alpha=1e-1;
beta=2;
k=0;
lambda=alpha^2*(n+k)-n;

% UT weights
Wm=zeros(2*n+1,1);
Wc=Wm;
Wm(1)=lambda/(lambda+n);
Wm(2:end)=1/(2*(lambda+n));
Wc(1)=Wm(1)+(1-alpha^2+beta);
Wc(2:end)=Wm(2:end);

% UKF operation
[xx_ukf_M,Px_ukf_M] = UKF(xs_M,1e4*Ps,ets,meas_S_M,tt_meas_S_M,str_S,R_S,Wm,Wc,lambda,opt);
fprintf('Final displacement: %.4e km, %.4e km/s\n',1*norm(xx_ukf_M(end,1:3)-xx_ref_M(end,1:3)),1*norm(xx_ukf_M(end,4:6)-xx_ref_M(end,4:6)))

% covariance position and velocity 3sigma
M3sr=zeros(1,length(tt_meas_S_M));
M3sv=zeros(1,length(tt_meas_S_M));

for i=1:length(tt_meas_S_M)
    M3sr(i)=3*sqrt(trace(Px_ukf_M(1:3,1:3,i)));
    M3sv(i)=3*sqrt(trace(Px_ukf_M(4:6,4:6,i)));
end

% 3sigma and displacements wrt reference SGP4 states
figure
subplot(2,1,1)
plot(5*(0:(length(tt_meas_S_M)-1))/60,sqrt(sum((xx_ukf_M(:,1:3)-xx_ref_M(:,1:3)).^2,2)),'x')
hold on
plot(5*(0:(length(tt_meas_S_M)-1))/60,M3sr,'r--','LineWidth',1)
grid minor
xlim([0 5*(length(tt_meas_S_M)-1)/60])
ylabel('$|\delta r|\, [km]$')
legend('$|\delta r|$','$3\sigma_r$','Orientation','horizontal','fontsize',25,'location','northeast')
title('Mango, ECI')
subplot(2,1,2)
plot(5*(0:(length(tt_meas_S_M)-1))/60,sqrt(sum((xx_ukf_M(:,4:6)-xx_ref_M(:,4:6)).^2,2)),'x')
hold on
plot(5*(0:(length(tt_meas_S_M)-1))/60,M3sv,'r--','LineWidth',1)
grid minor
xlim([0 5*(length(tt_meas_S_M)-1)/60])
ylabel('$|\delta v|\, [km/s]$')
xlabel(sprintf('Minutes past %s',cspice_et2utc(tt_meas_S_M(1),'C',0)))
legend('$|\delta v|$','$3\sigma_v$','Orientation','horizontal','fontsize',25,'location','northeast')


%% EX2
close all; clc;

% Tango TLE and satrec
T_longstr1='1 36827U 10028F   10224.22753605  .00278492  00000-0  82287-1 0  9996';
T_longstr2='2 36827 098.2797 049.5751 0044602 022.4408 337.8871 14.40890217    55';
T_rec=twoline2rv(T_longstr1, T_longstr2, typerun,'e', opsmode, whichconst);

[year,mon,day,hr,min,sec]=invjday(T_rec.jdsatepoch,T_rec.jdsatepochf);
et_T_TLE_str=sprintf('%d-%02d-%02d %02d:%02d:%02.6f',[year,mon,day,hr,min,sec]);

% Tango reference TLE epoch
et_T_TLE=cspice_str2et(et_T_TLE_str);

fprintf('Tango satellite ID: %d\n',T_rec.satnum);
fprintf('TLE reference epoch: UTC %s\n\n',et_T_TLE_str);

% nutation correction
arcsec2rad=pi/(180*3600);
ddpsi=-0.073296*arcsec2rad;
ddeps=-0.009373*arcsec2rad;

[~,rt_M,vt_M]=sgp4(M_rec,(et0-et_M_TLE)/60);
[~,rt_T,vt_T]=sgp4(T_rec,(et0-et_T_TLE)/60);
at=[0 0 0]';
ttt=cspice_unitim(et0,'ET','TDT')/cspice_jyear()/100;

% Mango eci state at et0
[r,v,~]=teme2eci(rt_M,vt_M,at,ttt,ddpsi,ddeps);
x0_M=[r', v'];

% Tango eci state at et0
[r,v,~]=teme2eci(rt_T,vt_T,at,ttt,ddpsi,ddeps);
x0_T=[r', v'];

% relative state at et0
x0_TM_eci=x0_T-x0_M;

% relative state at ets
xs_TM_eci=(xs_T-xs_M)';

% Mango mean motion
n_M=sqrt(mu/norm(x0_M(1:3))^3);

% state rotation matrix at et0 wrt Mango LVLH
RRM0=LVLH_rotm(x0_M,et0);

% relative LVLH state at et0
xr0_TM=(RRM0*x0_TM_eci')';

% CW propagate relative state
xr_TM=CW(tt,et0,xr0_TM,n_M);

% FFRF noise matrix
str_FFRF='FFRF';
R_FFRF=diag([deg2rad(1) deg2rad(1) 1e-5].^2);

% FFRF anticipated measurements
exp_meas_FFRF_0n=zeros(3,length(tt));
exp_meas_FFRF_0n(1,:)=atan2(xr_TM(:,2),xr_TM(:,1));
exp_meas_FFRF_0n(2,:)=asin(xr_TM(:,3)./sqrt(sum(xr_TM(:,1:3).^2,2)));
exp_meas_FFRF_0n(3,:)=sqrt(sum(xr_TM(:,1:3).^2,2));

% FFRF measurements
exp_meas_FFRF=(mvnrnd(exp_meas_FFRF_0n',R_FFRF))';

% FFRF measurement allocation
FFRF_id=find(tt==(VW_S_M(2,1)))+(1:(20*60/tm)+1);
tt_FFRF=tt(FFRF_id);
fprintf('%s measurement window from %s to %s\n',str_FFRF,cspice_et2utc(tt_FFRF(1),'C',1),cspice_et2utc(tt_FFRF(end),'C',1))
fprintf('Total %s measurements %.0f\n\n',str_FFRF,length(tt_FFRF))

% FFRF recorded measurements
meas_FFRF=exp_meas_FFRF(:,FFRF_id);

% Tango apparent motion wrt Mango
figure
plot(rad2deg(meas_FFRF(1,:)),rad2deg(meas_FFRF(2,:)),'rx')
hold on
plot(rad2deg(exp_meas_FFRF_0n(1,FFRF_id)),rad2deg(exp_meas_FFRF_0n(2,FFRF_id)),'bo')
xlim([-125 -65])
ylim([-15 -2])
grid minor
xlabel('$Az\,[deg]$')
ylabel('$El\,[deg]$')
title('Tango - Mango')

% FFRF measurements
figure
subplot(3,1,1)
plot(5*(0:(length(tt_FFRF)-1))/60,rad2deg(meas_FFRF(1,:)),'rx')
grid minor
xlim([1 5*(length(tt_FFRF)-1)/60])
ylabel('$Az\,[deg]$')
title(str_FFRF)
subplot(3,1,2)
plot(5*(0:(length(tt_FFRF)-1))/60,rad2deg(meas_FFRF(2,:)),'rx')
grid minor
xlim([1 5*(length(tt_FFRF)-1)/60])
ylabel('$El\,[deg]$')
subplot(3,1,3)
plot(5*(0:(length(tt_FFRF)-1))/60,(meas_FFRF(3,:)),'rx')
grid minor
xlim([1 5*(length(tt_FFRF)-1)/60])
ylabel('$r\,[km]$')
% xlabel('Measurements\,[-]')
xlabel(sprintf('Minutes past %s',cspice_et2utc(tt_FFRF(1),'C',0)))

% recorded reference state
xr_ref=xr_TM(FFRF_id,:);

% state rotation matrix at ets wrt Mango LVLH
RRMs=LVLH_rotm(xs_M,ets);

% LVLH relative UKF initial conditions
xrs_TM=RRMs*xs_TM_eci'; % relative mean state
Prs=2*RRMs*Ps*RRMs';    % relative covariance

% UKF computation
[xr_ukf,Pr_ukf]=UKF_FFRF(xrs_TM,Prs,ets,meas_FFRF,tt_FFRF,R_FFRF,Wm,Wc,lambda,n_M);

fprintf('Final displacement: %.4e km, %.4e km/s\n',1*norm(xr_ukf(end,1:3)-xr_ref(end,1:3)),1*norm(xr_ukf(end,4:6)-xr_ref(end,4:6)))

% relative covariance 3sigma allocation
r3sr=zeros(1,length(tt_FFRF));
r3sv=zeros(1,length(tt_FFRF));

% relative covariance 3sigma computation
for i=1:length(tt_FFRF)
    r3sr(i)=3*sqrt(trace(Pr_ukf(1:3,1:3,i)));
    r3sv(i)=3*sqrt(trace(Pr_ukf(4:6,4:6,i)));
end

% 3sigma and displacements wrt reference CW propagated states
figure
subplot(2,1,1)
plot(5*(0:(length(tt_FFRF)-1))/60,sqrt(sum((xr_ukf(:,1:3)-xr_ref(:,1:3)).^2,2)),'x')
hold on
plot(5*(0:(length(tt_FFRF)-1))/60,r3sr,'r--','LineWidth',1)
grid minor
xlim([1 5*(length(tt_FFRF)-1)/60])
ylabel('$|\delta r|\, [km]$')
legend('$|\delta r|$','$3\sigma_r$','Orientation','horizontal','fontsize',25,'location','northeast')
title('Tango-Mango, LVLH')
subplot(2,1,2)
plot(5*(0:(length(tt_FFRF)-1))/60,sqrt(sum((xr_ukf(:,4:6)-xr_ref(:,4:6)).^2,2)),'x')
hold on
plot(5*(0:(length(tt_FFRF)-1))/60,r3sv,'r--','LineWidth',1)
grid minor
xlim([1 5*(length(tt_FFRF)-1)/60])
ylabel('$|\delta v|\, [km/s]$')
legend('$|\delta v|$','$3\sigma_v$','Orientation','horizontal','fontsize',25,'location','northeast')
% xlabel('Measurements\, [-]')
xlabel(sprintf('Minutes past %s',cspice_et2utc(tt_FFRF(1),'C',0)))

%% EX3
close all; clc;

% Mango covariance propagation window
tt_T=[tt_meas_S_M(end) tt_FFRF];

XC=sqrtm((n+lambda)*Px_ukf_M(:,:,end));

% propagated sigma points allocation
YM_FFRF=zeros(length(tt_FFRF),n,2*n+1);
for i=1:2*n+1

    % sigma points computation
    if i==1
        XM=xx_ukf_M(end,:)';
    elseif i<=n+1
        XM=xx_ukf_M(end,:)'+XC(:,i-1);
    else
        XM=xx_ukf_M(end,:)'-XC(:,i-1-n);
    end    

    % sigma points propagation
    [~,Y]=ode113(@(t,x) KEPJ2(t,x,mu),tt_T,XM,opt);
    YM_FFRF(:,:,i)=Y(2:end,:);
end

% Mango UT mean states
xx_M_FFRF=sum(YM_FFRF.*reshape(Wm,1,1,[]),3);

% eci covariances allocations
Px_M_FFRF=zeros(n,n,length(tt_FFRF));   % absolute Mango covariance
Px_T_FFRF=Px_M_FFRF;                    % absolute Tango covariance
Pr_TM_FFRF=Px_M_FFRF;                   % relative covariance

% 3sigma allocation
T3sr_FFRF=zeros(1,length(tt_FFRF)); % Tango
T3sv_FFRF=T3sr_FFRF;
M3sr_FFRF=T3sr_FFRF;                % Mango
M3sv_FFRF=T3sv_FFRF;
r3sr_FFRF=T3sr_FFRF;                % relative Tango-Mango
r3sv_FFRF=T3sv_FFRF;

for i=1:length(tt_FFRF)

    % absolute Mango UT covariances
    for j=1:2*n+1
        Px_M_FFRF(:,:,i)=Px_M_FFRF(:,:,i)+Wc(j)*(xx_M_FFRF(i,:)-YM_FFRF(i,:,j))'*(xx_M_FFRF(i,:)-YM_FFRF(i,:,j));
    end

    % eci to lvlh stare rotation matrix
    R_eci2lvlh=LVLH_rotm(xx_M_FFRF(i,:),tt_FFRF(i));
    
    % relative eci Tango-Mango covariance
    Pr_TM_FFRF(:,:,i)=R_eci2lvlh\Pr_ukf(:,:,i)/R_eci2lvlh';

    % absolute Tango covariance
    Px_T_FFRF(:,:,i)=Px_M_FFRF(:,:,i)+Pr_TM_FFRF(:,:,i);

    % Tango 3sigma computation
    T3sr_FFRF(i)=3*sqrt(trace(Px_T_FFRF(1:3,1:3,i)));
    T3sv_FFRF(i)=3*sqrt(trace(Px_T_FFRF(4:6,4:6,i)));

    % Mango 3sigma computation
    M3sr_FFRF(i)=3*sqrt(trace(Px_M_FFRF(1:3,1:3,i)));
    M3sv_FFRF(i)=3*sqrt(trace(Px_M_FFRF(4:6,4:6,i)));

    % relative Tango-Mango 3sigma computation
    r3sr_FFRF(i)=3*sqrt(trace(Pr_TM_FFRF(1:3,1:3,i)));
    r3sv_FFRF(i)=3*sqrt(trace(Pr_TM_FFRF(4:6,4:6,i)));
end

% absolute Mango 3sigma and relative 3sigma
figure
subplot(2,1,1)
semilogy(5*(0:length(tt_FFRF)-1)/60,M3sr_FFRF,'-')
hold on
semilogy(5*(0:length(tt_FFRF)-1)/60,r3sr_FFRF,'-')
grid minor
xlim([1 5*(length(tt_FFRF)-1)/60])
ylabel('$3\sigma_r\, [km]$')
legend('$3\sigma_{r,M}$','$3\sigma_{r,TM}$','location','east','fontsize',25)
title('Tango, ECI')
subplot(2,1,2)
semilogy(5*(0:length(tt_FFRF)-1)/60,M3sv_FFRF,'-')
hold on
semilogy(5*(0:length(tt_FFRF)-1)/60,r3sv_FFRF,'-')
grid minor
xlim([1 5*(length(tt_FFRF)-1)/60])
ylabel('$3\sigma_v\, [km/s]$')
legend('$3\sigma_{v,M}$','$3\sigma_{v,TM}$','location','east','fontsize',25)
xlabel(sprintf('Minutes past %s',cspice_et2utc(tt_FFRF(1),'C',0)))

% absoluta Tango 3sigma
figure
subplot(2,1,1)
plot(5*(0:length(tt_FFRF)-1)/60,T3sr_FFRF,'-')
grid minor
xlim([1 5*(length(tt_FFRF)-1)/60])
ylabel('$3\sigma_r\, [km]$')
legend('$3\sigma_{r,T}$','location','northwest','fontsize',25)
title('Tango, ECI')
subplot(2,1,2)
plot(5*(0:length(tt_FFRF)-1)/60,T3sv_FFRF,'-')
grid minor
xlim([1 5*(length(tt_FFRF)-1)/60])
ylabel('$3\sigma_v\, [km/s]$')
legend('$3\sigma_{v,T}$','location','northwest','fontsize',25)
xlabel(sprintf('Minutes past %s',cspice_et2utc(tt_FFRF(1),'C',0)))

%% User plot setting removal
set(0,'DefaultAxesFontSize','remove');
set(0,'DefaultLegendFontSize','remove');
set(0,'DefaultLineLineWidth','remove');
set(groot,'defaultTextInterpreter','remove');
set(groot,'defaultLegendInterpreter','remove');
cspice_kclear;


%% Functions

function [rhs] = KEP(~,x,mu)
%     Computes keplerian equations of motion RHS
%     Example: [rhs] = KEP(t, x, mu)
%     INPUTS:
%         t  [1x1]  reference epoch (omissible)
%         x  [1x6]  state
%         mu [1x1]  KEP gravitational constant
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
    rhs(1:3)=[vx vy vz].';
    rhs(4:6)=-mu.*rr./norm(rr)^3;
end

function [rhs] = KEPJ2(t,x,mu)
%     Computes keplerian equations of motion RHS with J2 disturbance
%     Example: [rhs] = KEPJ2(t, x, mu)
%     INPUTS:
%         t  [1x1]  reference epoch
%         x  [1x6]  state
%         mu [1x1]  KEP gravitational constant
%     OUTPUTS:
%         rhs [6x1] equations of motion and STM propagation RHS

    % State verticality check
    [m,n]=size(x);
    if n>m
        x=x';
    end
    
    % ECI to ECEF rotation matrix
    eci2ecef=cspice_pxform('J2000','ITRF93',t);

    % Earth radius
    radii=cspice_bodvrd('Earth','RADII',3);
    Re=radii(1);

    % J2 coefficient
    j2=0.0010826269;

    % State Variables Extraction
    rx=x(1);
    ry=x(2);
    rz=x(3);
    vx=x(4);
    vy=x(5);
    vz=x(6);

    rr=[rx ry rz].';

    % ECEF state
    rr_ecef=eci2ecef*rr;

    % J2 ECEF acceleration
    aj2=3/2*mu*j2*rr_ecef/(norm(rr_ecef)^3)*((Re/norm(rr_ecef))^2).*(5*(rr_ecef(3)/norm(rr_ecef))^2-[1;1;3]);

    % J2 ECI acceleration
    aj2eci=eci2ecef'*aj2;
    
    % RHS assembly
    rhs=zeros(6,1);
    rhs(1:3)=[vx vy vz].';
    rhs(4:6)=-mu.*rr./norm(rr)^3+aj2eci;
end

function [xx] = CW(tt,t0,x0,nm)
%     Computes Clohessy-Wiltshire equations of motion RHS
%     Example: [rhs] = CW(tt, t0, x0, nm)
%     INPUTS:
%         tt [1xn]  integration epochs
%         t0 [1x1]  initial reference epoch
%         x0 [6x1]  relative state
%         nm [1x1]  mean angular velocyty 
%     OUTPUTS:
%         xx [nx6]  propagated state

    % State verticality check
    [m,n]=size(x0);
    if n>m
        x0=x0';
    end

    % STM submatrices
    Phi_rr=@(t) [4-3*cos(nm*(t-t0)) 0 0;
                 6*(sin(nm*(t-t0))-nm*(t-t0)) 1 0;
                 0 0 cos(nm*(t-t0))];

    Phi_rv=@(t) [sin(nm*(t-t0))/nm 2*(1-cos(nm*(t-t0)))/nm 0;
                 2*(cos(nm*(t-t0))-1)/nm (4*sin(nm*(t-t0))-3*nm*(t-t0))/nm 0;
                 0 0 sin(nm*(t-t0))/nm];

    Phi_vr=@(t) [3*nm*sin(nm*(t-t0)) 0 0;
                 6*nm*(cos(nm*(t-t0))-1) 0 0;
                 0 0 -nm*sin(nm*(t-t0))];

    Phi_vv=@(t) [cos(nm*(t-t0)) 2*sin(nm*(t-t0)) 0;
                 -2*sin(nm*(t-t0)) 4*cos(nm*(t-t0))-3 0;
                 0 0 cos(nm*(t-t0))];

    % propagated states allocation
    xx=zeros(length(tt),length(x0));

    for i=1:length(tt)
        % propagated states computation
        r=Phi_rr(tt(i))*x0(1:3)+Phi_rv(tt(i))*x0(4:6);
        v=Phi_vr(tt(i))*x0(1:3)+Phi_vv(tt(i))*x0(4:6);
        xx(i,:)=[r' v'];
    end
end

function [state,meas_n,tt_meas,vis_id]=station_measurements(satrec,station_name,tt,et_TLE)
    %     Computes Az,El,r measurements with respect to selected stations
    %     Example: [meas_n, tt_meas, vis_id] = station_measurements(satrec, station_name, tt, et_TLE)
    %     INPUTS:
    %         satrec       [struct]  satellite satrec
    %         station_name [str]  name of selected station
    %         tt           [1xn]  measurements acquisition epochs
    %         et_TLE       [1x1]  reference TLE epoch
    %     OUTPUTS:
    %         state        [nx6] unperturbed measurements associated states
    %         meas_n       [3xn] [Az,El,r] measurements [rad,rad,km]
    %         tt_meas      [1xm] epochs of filtered measurements
    %         vis_id       [1xn] unfiltered measurements identifier

    % sgp4 time grid
    tt_sgp4=(tt-et_TLE)/60;
    arcsec2rad=pi/(180*3600);

    % nutation correction
    ddpsi=-0.073296*arcsec2rad;
    ddeps=-0.009373*arcsec2rad;

    % measurements allocation
    meas=zeros(3,length(tt));

    % station selection
    if isequal(station_name,'KOUROU')
        R=diag([deg2rad(100*1e-3) deg2rad(100*1e-3) 0.01].^2);
        min_el=(10);
    elseif isequal(station_name,'SVALBARD')
        R=diag([deg2rad(125*1e-3) deg2rad(125*1e-3) 0.01].^2);
        min_el=(5);
    else
        error('Incorrect station name')
    end

    % unperturbed states allocation
    state=zeros(length(tt),6);

    for i=1:length(tt)

        % sgp4 state computation 
        [~,rt,vt]=sgp4(satrec,tt_sgp4(i));
        at=[0 0 0]';
        ttt=cspice_unitim(tt(i),'ET','TDT')/cspice_jyear()/100;
        [r,v,~]=teme2eci(rt,vt,at,ttt,ddpsi,ddeps);
        xx_sgp4=[r', v'];

        % spg4 computed unperturbed states
        state(i,:)=xx_sgp4;
        
        % perfect measurement acquisition
        [meas(1,i),meas(2,i),meas(3,i)]=antenna_pointing(station_name,tt(i),xx_sgp4);

    end

    % noise disturbance addition
    meas_n=meas+mvnrnd([0 0 0],R,length(tt))';

    % filtering of unfit measurements
    vis_id=meas_n(2,:)>deg2rad(min_el);
    state=state(vis_id,:);
    tt_meas=tt(vis_id);
    meas_n=meas_n(:,vis_id);

    fprintf('Total %s measurements %.0f, %.0f measurements filtered out\n',station_name,length(tt_meas),length(tt)-length(tt_meas))
end

function [az,el,range] = antenna_pointing(station_name, tt, xx)
    %     Computes Az,El,r conversion of x(t) with respect to selected stations
    %     Example: [az, el, range] = antenna_pointing(station_name, tt, xx)
    %     INPUTS:
    %         station_name [str]  name of selected station
    %         tt           [1xn]  reference epochs
    %         xx           [nx6]  reference states
    %     OUTPUTS:
    %         az      [1x1] azimuth [rad]
    %         el      [1x1] elevation [rad]
    %         r       [1x1] range [km]

    % horizontality check
    if size(xx,1)~=length(tt)
        xx=xx';
    end

    % station selection
    if ~isequal(station_name,'KOUROU') && ~isequal(station_name,'SVALBARD')
        error('Incorrect station name')
    end
    station_topo = [station_name, '_TOPO'];

    % conversion allocation
    az=zeros(1,length(tt));
    el=az;
    range=az;

    for i=1:length(tt)

        % ECI to ECEF TOPO rotation
        eci2topo=cspice_pxform('J2000', station_topo, tt(i));

        % ECEF TOPO station position
        rr_station_topo=cspice_spkpos(station_name, tt(i), 'J2000', 'NONE', 'EARTH');

        % relative ECEF TOPO reference position
        xx_topo = eci2topo*(xx(i,1:3)' - rr_station_topo);
        
        % az,el,r computation
        az(i)=atan2(xx_topo(2),xx_topo(1));
        el(i)=asin(xx_topo(3)/norm(xx_topo));
        range(i)=norm(xx_topo);
    end
end

function [xx,Px] = UKF(x0,P0,t0,meas,tt_meas,station_name,R,Wm,Wc,lambda,options)
    %     Computes UKF state and covariance estimates for KEPJ2 dynamics
    %     Example: [xx,Px] = UKF(x0, P0, t0, meas, tt_meas, station_name, R, Wm, Wc, lambda, options)
    %     INPUTS:
    %         x0           [1x6] initialization state
    %         P0           [6x6] initialization covariance
    %         t0           [1x1] initialization epoch
    %         meas         [3xn] measurements
    %         tt_meas      [1xn] measurements associated epochs
    %         station_name [str] station name
    %         R            [3x3] measurements associated noise matrix
    %         Wm           [1x13] UT mean weights
    %         Wc           [1x13] UT covariance weights
    %         lambda       [1x1]  UT lambda parameter
    %         options      [struct] integrator settings
    %     OUTPUTS:
    %         xx           [nx6] filtered states
    %         Px           [6x6xn] filtered covariances
    
    % number of state elements
    n=length(x0);

    % Earth gravitational parameter
    mu=cspice_bodvrd('Earth','GM',1);

    % filtered states and covariances initialization
    xx=zeros(length(tt_meas)+1,n);
    Px=zeros(n,n,length(tt_meas)+1);

    % propagation time grid
    tt=[t0 tt_meas];

    % initialization
    xx(1,:)=x0;
    Px(:,:,1)=P0;

    for i=1:length(tt_meas)
        
        % previous a priori state estimate
        x_prev_post=xx(i,:)';

        % allocations
        X=zeros(2*n+1,n);
        Y=zeros(size(meas,1),2*n+1);
        XC=sqrtm((n+lambda)*Px(:,:,i));

        % UT state and measurements propagation
        for j=1:2*n+1
            if j==1
                X_prev=x_prev_post;
            elseif j<=n+1
                X_prev=x_prev_post+XC(:,j-1);
            else
                X_prev=x_prev_post-XC(:,j-1-n);
            end   
            [~,XX]=ode113(@(t,x) KEPJ2(t,x,mu),[tt(i) tt(i+1)],X_prev,options);
            X(j,:)=XX(end,:);
            [Y(1,j),Y(2,j),Y(3,j)]=antenna_pointing(station_name,tt_meas(i),X(j,:));
        end
        
        % a priori mean state and mean measurement
        x_prior=sum(Wm.*X);
        y=sum(Wm'.*Y,2);

        % allocations
        P_prior=zeros(length(x_prior),length(x_prior));
        Pyy=zeros(size(meas,1),size(meas,1));
        Pxy=zeros(length(x_prior),size(meas,1));

        % a priori covariance, measurements covariance, meas-state cross
        % covariance
        for j=2*n+1:-1:1
            P_prior=P_prior+Wc(j).*(X(j,:)-x_prior)'*(X(j,:)-x_prior);
            Pyy=Pyy+Wc(j).*(Y(:,j)-y)*(Y(:,j)-y)';
            Pxy=Pxy+Wc(j).*(X(j,:)-x_prior)'*(Y(:,j)-y)';
        end
        Pyy=Pyy+R;

        % kalman gain
        K=Pxy/Pyy;

        % a posteriori correction
        xx(i+1,:)=x_prior+(K*(meas(:,i)-y))';
        P=P_prior-K*Pyy*K';

        % symmetry check
        Px(:,:,i+1)=0.5*(P+P');
    end

    xx=xx(2:end,:);
    Px=Px(:,:,2:end);
end

function [xx,Px] = UKF_FFRF(x0,P0,t0,meas,tt_meas,R,Wm,Wc,lambda,nm)
    %     Computes UKF state and covariance estimates for CW dymanics
    %     Example: [xx,Px] = UKF(x0, P0, t0, meas, tt_meas,, R, Wm, Wc, lambda, options)
    %     INPUTS:
    %         x0           [1x6] initialization state
    %         P0           [6x6] initialization covariance
    %         t0           [1x1] initialization epoch
    %         meas         [3xn] measurements
    %         tt_meas      [1xn] measurements associated epochs
    %         R            [3x3] measurements associated noise matrix
    %         Wm           [1x13] UT mean weights
    %         Wc           [1x13] UT covariance weights
    %         lambda       [1x1]  UT lambda parameter
    %         nm           [1x1] mean angular motion
    %     OUTPUTS:
    %         xx           [nx6] filtered states
    %         Px           [6x6xn] filtered covariances
    
    % number of state elements
    n=length(x0);

    % filtered states and covariances initialization
    xx=zeros(length(tt_meas)+1,n);
    Px=zeros(n,n,length(tt_meas)+1);

    % propagation time grid
    tt=[t0 tt_meas];

    % initialization
    xx(1,:)=x0;
    Px(:,:,1)=P0;


    for i=1:length(tt_meas)
        
        % previous a priori state estimate
        x_prev_post=xx(i,:)';

        % allocations
        X=zeros(2*n+1,n);
        Y=zeros(size(meas,1),2*n+1);
        XC=sqrtm((n+lambda)*Px(:,:,i));

        % UT state and measurements propagation
        for j=1:2*n+1
            if j==1
                X_prev=x_prev_post;
            elseif j<=n+1
                X_prev=x_prev_post+XC(:,j-1);
            else
                X_prev=x_prev_post-XC(:,j-1-n);
            end   
            XX=CW(tt(i+1),tt(i),X_prev,nm);
            X(j,:)=XX;
            Y(1,j)=atan2(XX(:,2),XX(:,1));
            Y(2,j)=asin(XX(:,3)/norm(XX(1:3)));
            Y(3,j)=norm(XX(1:3));
        end

        % a priori mean state and mean measurement
        x_prior=sum(Wm.*X);
        y=sum(Wm'.*Y,2);

        % allocations
        P_prior=zeros(length(x_prior),length(x_prior));
        Pyy=zeros(size(meas,1),size(meas,1));
        Pxy=zeros(length(x_prior),size(meas,1));

        % a priori covariance, measurements covariance, meas-state cross
        % covariance
        for j=2*n+1:-1:1
            P_prior=P_prior+Wc(j).*(X(j,:)-x_prior)'*(X(j,:)-x_prior);
            Pyy=Pyy+Wc(j).*(Y(:,j)-y)*(Y(:,j)-y)';
            Pxy=Pxy+Wc(j).*(X(j,:)-x_prior)'*(Y(:,j)-y)';
        end
        Pyy=Pyy+R;

        % kalman gain
        K=Pxy/Pyy;

        % a posteriori correction
        xx(i+1,:)=x_prior+(K*(meas(:,i)-y))';
        P=P_prior-K*Pyy*K';

        % symmetry check
        Px(:,:,i+1)=0.5*(P+P');
    end

    xx=xx(2:end,:);
    Px=Px(:,:,2:end);
end

function RR = LVLH_rotm(x,t)
    %     Computes state rotation matrix associated to x according to KEPJ2
    %     dynamics
    %     Example: RR = LVLH_rotm(x, t)
    %     INPUTS:
    %         x  [6x1]  reference state
    %         t  [1x1]  reference epoch
    %     OUTPUTS:
    %         RR [6x6] state rotation matrix

    % verticality check
    [m,n]=size(x);
    if n>m
        x=x';
    end
    r=x(1:3);
    v=x(4:6);

    % lvlh unit vectors
    i=r/norm(r);
    k=cross(r,v)/norm(cross(r,v));
    j=cross(k,i);

    % rotation matrix
    R=[i'; j'; k'];

    % Earth gravitational parameter
    mu=cspice_bodvrd('Earth','GM',1);

    % J2 acceleration computation
    rhs=KEPJ2(t,x,mu)-KEP(t,x,mu);
    aj=rhs(4:6);

    % lvlh unit vectors derivatives
    idot=(v-dot(v,i)*i)/norm(r);
    kdot=(cross(r,aj)-dot(cross(r,aj),k)*k)/norm(cross(r,v));
    jdot=cross(kdot,i)+cross(k,idot);

    % rotation matrix derivative
    Rdot=[idot'; jdot'; kdot'];

    % state rotation matrix
    RR=[R zeros(3,3);
        Rdot R];
end