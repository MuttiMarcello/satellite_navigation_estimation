% Spacecraft Guidance And Navigation
% Academic Year 2023/2024
% Assignment #2 exercise #2
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
etf=cspice_str2et('2010-08-12 11:00:00.000 UTC');

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

% Integrator settings
AbsTol=1e-12;   % Absolute tolerance
RelTol=1e-12;   % Relative tolerance
opt=odeset('AbsTol',AbsTol,'RelTol',RelTol);

% Mango state at et0
[~,xx0]=ode113(@(t,x) KEP(t,x,mu),[ets et0],xs_M,opt);
x0_M=xx0(end,:)';

% Time grid
tm=60;  % sampling time
tt=et0:tm:etf;

% keplerian Mango propagated states
[~,xx_M]=ode113(@(t,x) KEP(t,x,mu),tt,x0_M,opt);

% kourou 
str_K='KOUROU';
min_el_K=10; % [deg] minimum elevation requirement
% anticipated keplerian measurements
[az_M_K, el_M_K, r_M_K] = antenna_pointing(str_K, tt, xx_M);

% svalbard
str_S='SVALBARD';
min_el_S=5;  % [deg] minimum elevation requirement
% anticipated keplerian measurements
[az_M_S, el_M_S, r_M_S] = antenna_pointing(str_S, tt, xx_M);

% kourou and svalbard Mango expected measurements
figure
subplot(3,2,1)
plot(1:length(tt),rad2deg(az_M_K),'o')
grid minor
xlim([1 length(tt)])
ylim([-180 180])
ylabel('$Az\,[deg]$')
title(str_K)
subplot(3,2,2)
plot(1:length(tt),rad2deg(az_M_S),'o')
grid minor
xlim([1 length(tt)])
ylim([-180 180])
title(str_S)
subplot(3,2,3)
plot(1:length(tt),rad2deg(el_M_K),'o')
grid minor
xlim([1 length(tt)])
ylim([-90 90])
ylabel('$El\,[deg]$')
subplot(3,2,4)
plot(1:length(tt),rad2deg(el_M_S),'o')
grid minor
xlim([1 length(tt)])
ylim([-90 90])
subplot(3,2,5)
plot(1:length(tt),(r_M_K),'o')
grid minor
xlim([1 length(tt)])
xlabel(sprintf('Minutes past %s',cspice_et2utc(et0,'C',0)))
ylabel('r\,[km]')
subplot(3,2,6)
plot(1:length(tt),(r_M_S),'o')
grid minor
xlim([1 length(tt)])
xlabel(sprintf('Minutes past %s',cspice_et2utc(et0,'C',0)))

% visibility identifiers
K_vis_id_M=el_M_K>deg2rad(min_el_K);
S_vis_id_M=el_M_S>deg2rad(min_el_S);

% anticipated visible keplerian measurements
meas_kep_K_M=[az_M_K(K_vis_id_M==1); el_M_K(K_vis_id_M==1); r_M_K(K_vis_id_M==1)];
meas_kep_S_M=[az_M_S(S_vis_id_M==1); el_M_S(S_vis_id_M==1); r_M_S(S_vis_id_M==1)];

% anticipated visibility epochs
tt_K_M=tt(K_vis_id_M==1);
tt_S_M=tt(S_vis_id_M==1);

figure

% kourou apparent motion
subplot(1,2,1)
plot(rad2deg(az_M_K(K_vis_id_M)),rad2deg(el_M_K(K_vis_id_M)),'o')
hold on
yline(min_el_K,'-.','LineWidth',0.5)
grid minor
xlim([-180 180])
ylim([0 round(max(rad2deg([el_M_K(K_vis_id_M), el_M_S(S_vis_id_M)])),-1)+5*(round(max(rad2deg([el_M_K(K_vis_id_M),el_M_S(S_vis_id_M)])),-1)<max(rad2deg([el_M_K(K_vis_id_M),el_M_S(S_vis_id_M)])))])
xlabel('$Az\, [deg]$')
ylabel('$El\, [deg]$')
legend('Expected Keplerian Measurements','Location','south','fontsize',22)
title(str_K)

% allocated kourou visibility windows
VW_K_M=[];
VW_K_M(1,:)=tt(strfind(K_vis_id_M,[0 1])+1);  % first row opening of window
VW_K_M(2,:)=tt(strfind(K_vis_id_M,[1 0]));    % second row closing of window

fprintf('Mango\n')
for i=1:size(VW_K_M,2)
    fprintf('%s expected visibility window %.0f from %s to %s\n',str_K,i,cspice_et2utc(VW_K_M(1,i),'C',3),cspice_et2utc(VW_K_M(2,i),'C',3))
end
fprintf('Total %s expected visibility: %.2f min, %.0f measurements\n',str_K,sum(VW_K_M(2,:)-VW_K_M(1,:))/60,sum(K_vis_id_M))

% svalbard apparent motion
subplot(1,2,2)
plot(rad2deg(az_M_S(S_vis_id_M)),rad2deg(el_M_S(S_vis_id_M)),'o')
hold on
yline(min_el_S,'-.','LineWidth',0.5)
grid minor
xlim([-180 180])
ylim([0 round(max(rad2deg([el_M_K(K_vis_id_M), el_M_S(S_vis_id_M)])),-1)+5*(round(max(rad2deg([el_M_K(K_vis_id_M),el_M_S(S_vis_id_M)])),-1)<max(rad2deg([el_M_K(K_vis_id_M),el_M_S(S_vis_id_M)])))])
xlabel('$Az\, [deg]$')
ylabel('$El\, [deg]$')
legend('Expected Keplerian Measurements','Location','south','fontsize',22)
title(str_S)

% allocated svalbard visibility windows
VW_S_M=[];
VW_S_M(1,:)=tt(strfind(S_vis_id_M,[0 1])+1);  % first row opening of window
VW_S_M(2,:)=tt(strfind(S_vis_id_M,[1 0]));    % second row closing of window

for i=1:size(VW_S_M,2)
    fprintf('%s expected visibility window %.0f from %s to %s\n',str_S,i,cspice_et2utc(VW_S_M(1,i),'C',3),cspice_et2utc(VW_S_M(2,i),'C',3))
end
fprintf('Total %s expected visibility: %.2f min, %.0f measurements\n\n',str_S,sum(VW_S_M(2,:)-VW_S_M(1,:))/60,sum(S_vis_id_M))

% Earth plot
Re=cspice_bodvrd('Earth','RADII',3);
[X,Y,Z]=sphere(30);

% earth rotation velocity
w=7.2921159e-5;

% ECEF kourou position
rr_K=1.01*cspice_spkpos(str_K, tt_K_M(1), 'J2000', 'NONE', 'EARTH');

% ECEF svalbard position
rr_S=1.01*cspice_spkpos(str_S, tt_S_M(1), 'J2000', 'NONE', 'EARTH');

% Mango motion representation
figure
plot3(xx_M(K_vis_id_M,1)'.*cos(w*(tt_K_M-tt_K_M(1)))+xx_M(K_vis_id_M,2)'.*sin(w*(tt_K_M-tt_K_M(1))),-xx_M(K_vis_id_M,1)'.*sin(w*(tt_K_M-tt_K_M(1)))+xx_M(K_vis_id_M,2)'.*cos(w*(tt_K_M-tt_K_M(1))),xx_M(K_vis_id_M,3),'o')
hold on
plot3(xx_M(S_vis_id_M,1)'.*cos(w*(tt_S_M-tt_S_M(1)))+xx_M(S_vis_id_M,2)'.*sin(w*(tt_S_M-tt_S_M(1))),-xx_M(S_vis_id_M,1)'.*sin(w*(tt_S_M-tt_S_M(1)))+xx_M(S_vis_id_M,2)'.*cos(w*(tt_S_M-tt_S_M(1))),xx_M(S_vis_id_M,3),'o')
E=surf(Re(1)*X,Re(2)*Y,Re(3)*Z);
axis equal
view([115 45])
grid on
legend(str_K,str_S,'location','north')
set(E,'FaceColor',[0.7 0.7 0.7],'Edgecolor',[0.5 0.5 0.5])

%% EX2
close all; clc;

typerun    = 'u';  % user-provided inputs to SGP4 Matlab function
opsmode    = 'a';  % afspc approach ('air force space command')
whichconst =  72;  % WGS72 constants (radius, gravitational parameter)

% Mango TLE and satrec
M_longstr1 = '1 36599U 10028B   10224.22752732 -.00000576  00000-0 -16475-3 0  9998';
M_longstr2 = '2 36599 098.2803 049.5758 0043871 021.7908 338.5082 14.40871350  8293';
M_rec = twoline2rv(M_longstr1, M_longstr2, typerun,'e', opsmode, whichconst);

[year,mon,day,hr,min,sec] = invjday(M_rec.jdsatepoch, M_rec.jdsatepochf);
et_M_TLE_str=sprintf('%d-%02d-%02d %02d:%02d:%02.6f', [year,mon,day,hr,min,sec]);

% Mango reference TLE epoch
et_M_TLE=cspice_str2et(et_M_TLE_str);

fprintf('Mango satellite ID: %d\n', M_rec.satnum);
fprintf('TLE reference epoch: UTC %s\n\n', et_M_TLE_str);

% kourou Mango measurements
[meas_K_M,tt_meas_K_M,K_vis_id_M]=station_measurements(M_rec,str_K,tt_K_M,et_M_TLE);

% svalbard Mango measurements
[meas_S_M,tt_meas_S_M,S_vis_id_M]=station_measurements(M_rec,str_S,tt_S_M,et_M_TLE);

% stations apparent motion
figure
subplot(1,2,1)
plot(rad2deg(meas_kep_K_M(1,:)),rad2deg(meas_kep_K_M(2,:)),'o')
hold on
plot(rad2deg(meas_K_M(1,:)),rad2deg(meas_K_M(2,:)),'x')
yline(min_el_K,'-.','LineWidth',0.5)
grid minor
xlim([-180 180])
ylim([0 round(max(rad2deg([meas_kep_K_M(2,:), meas_K_M(2,:)])),-1)+5*(round(max(rad2deg([meas_kep_K_M(2,:), meas_K_M(2,:)])),-1)<max(rad2deg([meas_kep_K_M(2,:), meas_K_M(2,:)])))])
xlabel('$Az\, [deg]$')
ylabel('$El\, [deg]$')
legend('Expected Keplerian Measurements','Measurements','Location','south','fontsize',22)
title(str_K)
subplot(1,2,2)
plot(rad2deg(meas_kep_S_M(1,:)),rad2deg(meas_kep_S_M(2,:)),'o')
hold on
plot(rad2deg(meas_S_M(1,:)),rad2deg(meas_S_M(2,:)),'x')
yline(min_el_S,'-.','LineWidth',0.5)
grid minor
xlim([-180 180])
ylim([0 round(max(rad2deg([meas_kep_S_M(2,:), meas_S_M(2,:)])),-1)+5*(round(max(rad2deg([meas_kep_S_M(2,:), meas_S_M(2,:)])),-1)<max(rad2deg([meas_kep_S_M(2,:), meas_S_M(2,:)])))])
xlabel('$Az\, [deg]$')
ylabel('$El\, [deg]$')
legend('Expected Keplerian Measurements','Measurements','Location','south','fontsize',22)
title(str_S)

% stations measurements
figure
subplot(3,2,1)
plot(1:length(tt_meas_K_M),rad2deg(meas_K_M(1,:)),'rx')
grid minor
xlim([1 length(tt_meas_K_M)])
ylabel('$Az\,[deg]$')
title(str_K)
subplot(3,2,3)
plot(1:length(tt_meas_K_M),rad2deg(meas_K_M(2,:)),'rx')
grid minor
xlim([1 length(tt_meas_K_M)])
ylabel('$El\,[deg]$')
subplot(3,2,5)
plot(1:length(tt_meas_K_M),(meas_K_M(3,:)),'rx')
grid minor
xlim([1 length(tt_meas_K_M)])
ylabel('$r\,[km]$')
xlabel('Measurements\,[-]')
subplot(3,2,2)
plot(1:length(tt_meas_S_M),rad2deg(meas_S_M(1,:)),'rx')
grid minor
xlim([1 length(tt_meas_S_M)])
title(str_S)
subplot(3,2,4)
plot(1:length(tt_meas_S_M),rad2deg(meas_S_M(2,:)),'rx')
grid minor
xlim([1 length(tt_meas_S_M)])
subplot(3,2,6)
plot(1:length(tt_meas_S_M),(meas_S_M(3,:)),'rx')
grid minor
xlim([1 length(tt_meas_S_M)])
xlabel('Measurements\,[-]')

%% EX3
close all; clc;

% Mango teme state at et0
[~,rt,vt]=sgp4(M_rec,(et0-et_M_TLE)/60);
at=[0 0 0]';
ttt=cspice_unitim(et_M_TLE,'ET','TDT')/cspice_jyear()/100;
arcsec2rad=pi/(180*3600);

% nutation correction
ddpsi=-0.073296*arcsec2rad;
ddeps=-0.009373*arcsec2rad;

% Mango eci state at et0
[r,v,~]=teme2eci(rt,vt,at,ttt,ddpsi,ddeps);
x0guess_M=[r', v'];

% Least Squares algorithm settings
lsqopt=optimoptions('lsqnonlin','Display', 'iter','Algorithm','levenberg-marquardt');

% KOUROU meas, kep motion
tt_lsq_K_M=[et0, tt_meas_K_M];
[x0_K_M,resnorm_K_M,res_K_M,~,~,~,jacobian]=lsqnonlin(@(x0) cost_fun(x0,tt_lsq_K_M,meas_K_M,str_K,opt), x0guess_M, [], [], lsqopt);

% solution uncertainty covariance
Jac = full(jacobian);
P_lsq_K_M=resnorm_K_M/(length(res_K_M)-length(x0guess_M)).*inv(Jac.'*Jac);

% KOUROU and SVALBARD meas, kep motion
tt_lsq_S_M=[et0, tt_meas_S_M];
[x0_KS_M,resnorm_KS_M,res_KS_M,~,~,~,jacobian]=lsqnonlin(@(x0) [cost_fun(x0,tt_lsq_S_M,meas_S_M,str_S,opt); cost_fun(x0,tt_lsq_K_M,meas_K_M,str_K,opt)], x0guess_M, [], [], lsqopt);

% solution uncertainty covariance
Jac = full(jacobian);
P_lsq_KS_M=resnorm_KS_M/(length(res_KS_M)-length(x0guess_M)).*inv(Jac.'*Jac);

% KOUROU and SVALBARD meas, kep + J2 perturbation
[x0_KS_J2_M,resnorm_KS_J2_M,res_KS_J2_M,~,~,~,jacobian]=lsqnonlin(@(x0) [cost_fun(x0,tt_lsq_S_M,meas_S_M,str_S,opt,1); cost_fun(x0,tt_lsq_K_M,meas_K_M,str_K,opt,1)], x0guess_M, [], [], lsqopt);

% solution uncertainty covariance
Jac = full(jacobian);
P_lsq_KS_J2_M=resnorm_KS_J2_M/(length(res_KS_J2_M)-length(x0guess_M)).*inv(Jac.'*Jac);

% stations measures noise matrices
R_K=diag([deg2rad(100*1e-3) deg2rad(100*1e-3) 0.01].^2);
R_S=diag([deg2rad(125*1e-3) deg2rad(125*1e-3) 0.01].^2);

% kourou, kep residuals
res_K_M=reshape(res_K_M,[3,length(res_K_M)/3]);
res_K_M=sqrtm(R_K)*res_K_M;

figure
subplot(3,1,1)
plot(1:length(tt_lsq_K_M)-1,abs(rad2deg(res_K_M(1,:))),'o')
grid minor
xlim([1 length(tt_lsq_K_M)-1])
ylabel('$|\delta Az|\, [deg]$')
title(sprintf('%s, KEP',str_K))
subplot(3,1,2)
plot(1:length(tt_lsq_K_M)-1,abs(rad2deg(res_K_M(2,:))),'o')
grid minor
xlim([1 length(tt_lsq_K_M)-1])
ylabel('$|\delta El|\, [deg]$')
subplot(3,1,3)
plot(1:length(tt_lsq_K_M)-1,abs(res_K_M(3,:)),'o')
grid minor
xlim([1 length(tt_lsq_K_M)-1])
ylabel('$|\delta r|\, [km]$')
xlabel('Measurements\,[-]')

% kourou and svalbard, kep residuals
res_KS_M=reshape(res_KS_M,[3,length(res_KS_M)/3]);
res_KS_M(:,1:length(tt_meas_S_M))=sqrtm(R_S)*res_KS_M(:,1:length(tt_meas_S_M));
res_KS_M(:,length(tt_meas_S_M)+1:end)=sqrtm(R_K)*res_KS_M(:,length(tt_meas_S_M)+1:end);

figure
subplot(3,1,1)
plot(1:(length(tt_lsq_K_M)+length(tt_lsq_S_M))-2,abs(rad2deg(res_KS_M(1,:))),'o')
grid minor
xlim([1 (length(tt_lsq_K_M)+length(tt_lsq_S_M))-2])
ylabel('$|\delta Az|\, [deg]$')
title(sprintf('%s and %s, KEP',str_K,str_S))
subplot(3,1,2)
plot(1:(length(tt_lsq_K_M)+length(tt_lsq_S_M))-2,abs(rad2deg(res_KS_M(2,:))),'o')
grid minor
xlim([1 (length(tt_lsq_K_M)+length(tt_lsq_S_M))-2])
ylabel('$|\delta El|\, [deg]$')
subplot(3,1,3)
plot(1:(length(tt_lsq_K_M)+length(tt_lsq_S_M))-2,abs(res_KS_M(3,:)),'o')
grid minor
xlim([1 (length(tt_lsq_K_M)+length(tt_lsq_S_M))-2])
ylabel('$|\delta r|\, [km]$')
xlabel('Measurements\,[-]')

% kourou and svalbard, kep + J2 pert residuals
res_KS_J2_M=reshape(res_KS_J2_M,[3,length(res_KS_J2_M)/3]);
res_KS_J2_M(:,1:length(tt_meas_S_M))=sqrtm(R_S)*res_KS_J2_M(:,1:length(tt_meas_S_M));
res_KS_J2_M(:,length(tt_meas_S_M)+1:end)=sqrtm(R_K)*res_KS_J2_M(:,length(tt_meas_S_M)+1:end);

figure
subplot(3,1,1)
plot(1:(length(tt_lsq_K_M)+length(tt_lsq_S_M))-2,abs(rad2deg(res_KS_J2_M(1,:))),'o')
grid minor
xlim([1 (length(tt_lsq_K_M)+length(tt_lsq_S_M))-2])
ylabel('$|\delta Az|\, [deg]$')
title(sprintf('%s and %s, KEP with J2',str_K,str_S))
subplot(3,1,2)
plot(1:(length(tt_lsq_K_M)+length(tt_lsq_S_M))-2,abs(rad2deg(res_KS_J2_M(2,:))),'o')
grid minor
xlim([1 (length(tt_lsq_K_M)+length(tt_lsq_S_M))-2])
ylabel('$|\delta El|\, [deg]$')
subplot(3,1,3)
plot(1:(length(tt_lsq_K_M)+length(tt_lsq_S_M))-2,abs(res_KS_J2_M(3,:)),'o')
grid minor
xlim([1 (length(tt_lsq_K_M)+length(tt_lsq_S_M))-2])
ylabel('$|\delta r|\, [km]$')
xlabel('Measurements\,[-]')

fprintf('Mango initial state estimation using %s measurements:\n[%.4f %.4f %.4f %.4f %.4f %.4f]\n',str_K,x0_K_M(1),x0_K_M(2),x0_K_M(3),x0_K_M(4),x0_K_M(5),x0_K_M(6))
fprintf('Least Square Estimation covariance trace: %.4e\n',trace(P_lsq_K_M))
fprintf('Position STD deviation: %.4e Velocity STD deviation: %.4e\n',sqrt(trace(P_lsq_K_M(1:3,1:3))),sqrt(trace(P_lsq_K_M(4:6,4:6))))
fprintf('Absolute residual means: [%.4e %.4e %.4e]\n\n',mean(rad2deg(abs(res_K_M(1,:)))),mean(rad2deg(abs(res_K_M(2,:)))),mean((abs(res_K_M(3,:)))))

fprintf('Mango initial state estimation using %s and %s measurements:\n[%.4f %.4f %.4f %.4f %.4f %.4f]\n',str_K,str_S,x0_KS_M(1),x0_KS_M(2),x0_KS_M(3),x0_KS_M(4),x0_KS_M(5),x0_KS_M(6))
fprintf('Least Square Estimation covariance trace: %.4e\n',trace(P_lsq_KS_M))
fprintf('Position STD deviation: %.4e Velocity STD deviation: %.4e\n',sqrt(trace(P_lsq_KS_M(1:3,1:3))),sqrt(trace(P_lsq_KS_M(4:6,4:6))))
fprintf('Absolute residual means: [%.4e %.4e %.4e]\n\n',mean(rad2deg(abs(res_KS_M(1,:)))),mean(rad2deg(abs(res_KS_M(2,:)))),mean((abs(res_KS_M(3,:)))))

fprintf('Mango initial state estimation using %s and %s measurements and J2 perturbation:\n[%.4f %.4f %.4f %.4f %.4f %.4f]\n',str_K,str_S,x0_KS_J2_M(1),x0_KS_J2_M(2),x0_KS_J2_M(3),x0_KS_J2_M(4),x0_KS_J2_M(5),x0_KS_J2_M(6))
fprintf('Least Square Estimation covariance trace: %.4e\n',trace(P_lsq_KS_J2_M))
fprintf('Position STD deviation: %.4e Velocity STD deviation: %.4e\n',sqrt(trace(P_lsq_KS_J2_M(1:3,1:3))),sqrt(trace(P_lsq_KS_J2_M(4:6,4:6))))
fprintf('Absolute residual means: [%.4e %.4e %.4e]\n\n',mean(rad2deg(abs(res_KS_J2_M(1,:)))),mean(rad2deg(abs(res_KS_J2_M(2,:)))),mean((abs(res_KS_J2_M(3,:)))))

%% EX5
close all; clc;

% Tango state at et0
[~,xx0]=ode113(@(t,x) KEP(t,x,mu),[ets et0],xs_T,opt);
x0_T=xx0(end,:)';

% keplerian Tango propagated states
[~,xx_T]=ode113(@(t,x) KEP(t,x,mu),tt,x0_T,opt);

% anticipated kourou measurements
[az_T_K, el_T_K, r_T_K]=antenna_pointing(str_K,tt,xx_T);

% anticipated svalbard measurements
[az_T_S, el_T_S, r_T_S]=antenna_pointing(str_S,tt,xx_T);

% visibility identifiers
K_vis_id_T=el_T_K>deg2rad(min_el_K);
S_vis_id_T=el_T_S>deg2rad(min_el_S);

% anticipated visible keplerian measurements
meas_kep_K_T=[az_T_K(K_vis_id_T==1); el_T_K(K_vis_id_T==1); r_T_K(K_vis_id_T==1)];
meas_kep_S_T=[az_T_S(S_vis_id_T==1); el_T_S(S_vis_id_T==1); r_T_S(S_vis_id_T==1)];

% anticipated visibility epochs
tt_K_T=tt(K_vis_id_T==1);
tt_S_T=tt(S_vis_id_T==1);

% Tango motion representation
figure
plot3(xx_T(K_vis_id_T,1)'.*cos(w*(tt_K_T-tt_K_T(1)))+xx_T(K_vis_id_T,2)'.*sin(w*(tt_K_T-tt_K_T(1))),-xx_T(K_vis_id_T,1)'.*sin(w*(tt_K_T-tt_K_T(1)))+xx_T(K_vis_id_T,2)'.*cos(w*(tt_K_T-tt_K_T(1))),xx_T(K_vis_id_T,3),'o')
hold on
plot3(xx_T(S_vis_id_T,1)'.*cos(w*(tt_S_T-tt_S_T(1)))+xx_T(S_vis_id_T,2)'.*sin(w*(tt_S_T-tt_S_T(1))),-xx_T(S_vis_id_T,1)'.*sin(w*(tt_S_T-tt_S_T(1)))+xx_T(S_vis_id_T,2)'.*cos(w*(tt_S_T-tt_S_T(1))),xx_T(S_vis_id_T,3),'o')
E=surf(Re(1)*X,Re(2)*Y,Re(3)*Z);
axis equal
view([115 45])
grid on
legend(str_K,str_S,'location','north')
set(E,'FaceColor',[0.7 0.7 0.7],'Edgecolor',[0.5 0.5 0.5])

% allocated kourou visibility windows
VW_K_T=[];
VW_K_T(1,:)=tt(strfind(K_vis_id_T,[0 1])+1);  % first row opening of window
VW_K_T(2,:)=tt(strfind(K_vis_id_T,[1 0]));    % second row closing of window

fprintf('Tango\n')
for i=1:size(VW_K_T,2)
    fprintf('%s expected visibility window %.0f from %s to %s\n',str_K,i,cspice_et2utc(VW_K_T(1,i),'C',3),cspice_et2utc(VW_K_T(2,i),'C',3))
end
fprintf('Total %s expected visibility: %.2f min, %.0f measurements\n',str_K,sum(VW_K_T(2,:)-VW_K_T(1,:))/60,sum(K_vis_id_T))

% allocated svalbard visibility windows
VW_S_T=[];
VW_S_T(1,:)=tt(strfind(S_vis_id_T,[0 1])+1);  % first row opening of window
VW_S_T(2,:)=tt(strfind(S_vis_id_T,[1 0]));    % second row closing of window

for i=1:size(VW_S_T,2)
    fprintf('%s expected visibility window %.0f from %s to %s\n',str_S,i,cspice_et2utc(VW_S_T(1,i),'C',3),cspice_et2utc(VW_S_T(2,i),'C',3))
end
fprintf('Total %s expected visibility: %.2f min, %.0f measurements\n\n',str_S,sum(VW_S_T(2,:)-VW_S_T(1,:))/60,sum(S_vis_id_T))

% Tango TLE and satrec
T_longstr1 = '1 36827U 10028F   10224.22753605  .00278492  00000-0  82287-1 0  9996';
T_longstr2 = '2 36827 098.2797 049.5751 0044602 022.4408 337.8871 14.40890217    55';
T_rec = twoline2rv(T_longstr1, T_longstr2, typerun,'e', opsmode, whichconst);

[year,mon,day,hr,min,sec] = invjday(T_rec.jdsatepoch, T_rec.jdsatepochf);
et_T_TLE_str=sprintf('%d-%02d-%02d %02d:%02d:%02.6f', [year,mon,day,hr,min,sec]);

% Tango reference TLE epoch
et_T_TLE=cspice_str2et(et_T_TLE_str);

fprintf('Tango satellite ID: %d\n', T_rec.satnum);
fprintf('TLE reference epoch: UTC %s\n\n', et_T_TLE_str);

% kourou Tango measurements
[meas_K_T,tt_meas_K_T,K_vis_id_T]=station_measurements(T_rec,str_K,tt_K_T,et_T_TLE);

% svalbard Tango measurements
[meas_S_T,tt_meas_S_T,S_vis_id_T]=station_measurements(T_rec,str_S,tt_S_T,et_T_TLE);

% stations apparent motion
figure
subplot(1,2,1)
plot(rad2deg(meas_kep_K_T(1,:)),rad2deg(meas_kep_K_T(2,:)),'o')
hold on
plot(rad2deg(meas_K_T(1,:)),rad2deg(meas_K_T(2,:)),'x')
yline(min_el_K,'-.','LineWidth',0.5)
grid minor
xlim([-180 180])
ylim([0 round(max(rad2deg([meas_kep_K_T(2,:), meas_K_T(2,:)])),-1)+5*(round(max(rad2deg([meas_kep_K_T(2,:), meas_K_T(2,:)])),-1)<max(rad2deg([meas_kep_K_T(2,:), meas_K_T(2,:)])))])
xlabel('$Az\, [deg]$')
ylabel('$El\, [deg]$')
legend('Expected Keplerian Measurements','Measurements','Location','south','fontsize',22)
title(str_K)
subplot(1,2,2)
plot(rad2deg(meas_kep_S_T(1,:)),rad2deg(meas_kep_S_T(2,:)),'o')
hold on
plot(rad2deg(meas_S_T(1,:)),rad2deg(meas_S_T(2,:)),'x')
yline(min_el_S,'-.','LineWidth',0.5)
grid minor
xlim([-180 180])
ylim([0 round(max(rad2deg([meas_kep_S_T(2,:), meas_S_T(2,:)])),-1)+5*(round(max(rad2deg([meas_kep_S_T(2,:), meas_S_T(2,:)])),-1)<max(rad2deg([meas_kep_S_T(2,:), meas_S_T(2,:)])))])
xlabel('$Az\, [deg]$')
ylabel('$El\, [deg]$')
legend('Expected Keplerian Measurements','Measurements','Location','south','fontsize',22)
title(str_S)

% stations measurements
figure
subplot(3,2,1)
plot(1:length(tt_meas_K_T),rad2deg(meas_K_T(1,:)),'rx')
grid minor
xlim([1 length(tt_meas_K_T)])
ylabel('$Az\,[deg]$')
title(str_K)
subplot(3,2,3)
plot(1:length(tt_meas_K_T),rad2deg(meas_K_T(2,:)),'rx')
grid minor
xlim([1 length(tt_meas_K_T)])
ylabel('$El\,[deg]$')
subplot(3,2,5)
plot(1:length(tt_meas_K_T),(meas_K_T(3,:)),'rx')
grid minor
xlim([1 length(tt_meas_K_T)])
ylabel('$r\,[km]$')
xlabel('Measurements\,[-]')
subplot(3,2,2)
plot(1:length(tt_meas_S_T),rad2deg(meas_S_T(1,:)),'rx')
grid minor
xlim([1 length(tt_meas_S_T)])
title(str_S)
subplot(3,2,4)
plot(1:length(tt_meas_S_T),rad2deg(meas_S_T(2,:)),'rx')
grid minor
xlim([1 length(tt_meas_S_T)])
subplot(3,2,6)
plot(1:length(tt_meas_S_T),(meas_S_T(3,:)),'rx')
grid minor
xlim([1 length(tt_meas_S_T)])
xlabel('Measurements\,[-]')

% Tango teme state at et0
[~,rt,vt]=sgp4(T_rec,(et0-et_T_TLE)/60);
at=[0 0 0]';
ttt=cspice_unitim(et_M_TLE,'ET','TDT')/cspice_jyear()/100;

% Tango eci state at et0
[r,v,~]=teme2eci(rt,vt,at,ttt,ddpsi,ddeps);
x0guess_T=[r', v'];

% KOUROU and SVALBARD meas, kep + J2 pert
tt_lsq_K_T=[et0, tt_meas_K_T];
tt_lsq_S_T=[et0, tt_meas_S_T];
[x0_KS_J2_T,resnorm_KS_J2_T,res_KS_J2_T,~,~,~,jacobian]=lsqnonlin(@(x0) [cost_fun(x0,tt_lsq_S_T,meas_S_T,str_S,opt,1); cost_fun(x0,tt_lsq_K_T,meas_K_T,str_K,opt,1)], x0guess_T, [], [], lsqopt);

% solution uncertainty covariance
Jac = full(jacobian);
P_lsq_KS_J2_T=resnorm_KS_J2_T/(length(res_KS_J2_T)-length(x0guess_T)).*inv(Jac.'*Jac);

% residuals
res_KS_J2_T=reshape(res_KS_J2_T,[3,length(res_KS_J2_T)/3]);
res_KS_J2_T(:,1:length(tt_meas_S_T))=sqrtm(R_S)*res_KS_J2_T(:,1:length(tt_meas_S_T));
res_KS_J2_T(:,length(tt_meas_S_T)+1:end)=sqrtm(R_K)*res_KS_J2_T(:,length(tt_meas_S_T)+1:end);

figure
subplot(3,1,1)
plot(1:(length(tt_lsq_K_T)+length(tt_lsq_S_T))-2,abs(rad2deg(res_KS_J2_T(1,:))),'o')
grid minor
xlim([1 (length(tt_lsq_K_T)+length(tt_lsq_S_T))-2])
ylabel('$|\delta Az|\, [deg]$')
title(sprintf('%s and %s, KEP with J2',str_K,str_S))
subplot(3,1,2)
plot(1:(length(tt_lsq_K_T)+length(tt_lsq_S_T))-2,abs(rad2deg(res_KS_J2_T(2,:))),'o')
grid minor
xlim([1 (length(tt_lsq_K_T)+length(tt_lsq_S_T))-2])
ylabel('$|\delta El|\, [deg]$')
subplot(3,1,3)
plot(1:(length(tt_lsq_K_T)+length(tt_lsq_S_T))-2,abs(res_KS_J2_T(3,:)),'o')
grid minor
xlim([1 (length(tt_lsq_K_T)+length(tt_lsq_S_T))-2])
ylabel('$|\delta r|\, [km]$')
xlabel('Measurements\,[-]')

fprintf('\nTango initial state estimation using %s and %s measurements and J2 perturbation:\n[%.4f %.4f %.4f %.4f %.4f %.4f]\n',str_K,str_S,x0_KS_J2_T(1),x0_KS_J2_T(2),x0_KS_J2_T(3),x0_KS_J2_T(4),x0_KS_J2_T(5),x0_KS_J2_T(6))
fprintf('Least Square Estimation covariance trace: %.4e\n',trace(P_lsq_KS_J2_T))
fprintf('Position STD deviation: %.4e Velocity STD deviation: %.4e\n\n',sqrt(trace(P_lsq_KS_J2_T(1:3,1:3))),sqrt(trace(P_lsq_KS_J2_T(4:6,4:6))))


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

function [meas_n,tt_meas,vis_id]=station_measurements(satrec,station_name,tt,et_TLE)
    %     Computes Az,El,r measurements with respect to selected stations
    %     Example: [meas_n, tt_meas, vis_id] = station_measurements(satrec, station_name, tt, et_TLE)
    %     INPUTS:
    %         satrec       [struct]  satellite satrec
    %         station_name [str]  name of selected station
    %         tt           [1xn]  measurements acquisition epochs
    %         et_TLE       [1x1]  reference TLE epoch
    %     OUTPUTS:
    %         meas_n       [3xn] [Az,El,r] measurements [rad,rad,km]
    %         tt_meas      [1xm] epochs of filtered measurements
    %         vis_id       [1xn] unfiltered measurements identifier

    % sgp4 time grid
    tt_sgp4=(tt-et_TLE)/60;
    arcsec2rad=pi/(180*3600);

    % nutation correction
    ddpsi=-0.073296*arcsec2rad;
    ddeps=-0.009373*arcsec2rad;

    % measurement allocation
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

    for i=1:length(tt)

        % sgp4 state computation
        [~,rt,vt]=sgp4(satrec,tt_sgp4(i));     
        at=[0 0 0]';
        ttt=cspice_unitim(tt(i),'ET','TDT')/cspice_jyear()/100;
        [r,v,~]=teme2eci(rt,vt,at,ttt,ddpsi,ddeps);
        xx_sgp4=[r', v'];
        
        % perfect measurement acquisition
        [meas(1,i),meas(2,i),meas(3,i)]=antenna_pointing(station_name,tt(i),xx_sgp4);

    end

    % noise disturbance addition
    meas_n=meas+mvnrnd([0 0 0],R,length(tt))';
    
    % filtering of unfit measurements
    vis_id=meas_n(2,:)>deg2rad(min_el);
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

function [res] = cost_fun(x0,tt,meas,station_name,opt,j2)
    %     Computes residuals for batch filter WLSQ computation
    %     Example: [res] = cost_fun(x0guess, tt, meas, station_name, opt, j2)
    %     INPUTS:
    %         x0           [1x6]   initial propagation condition
    %         tt           [1xn+1] initial and measurements reference epochs
    %         meas         [3xn]   az,el,r measurements
    %         station_name [str]   name of selected station
    %         opt          [struct]integrator settings
    %         j2           [1x1]   j2 pert activation flag
    %     OUTPUTS:
    %         res          [3nx1] residuals [rad,km]

    % j2 pert null by default
    if nargin<6
        j2=0;
    end

    % station selection
    if isequal(station_name,'KOUROU')
        R=diag([deg2rad(100*1e-3) deg2rad(100*1e-3) 0.01].^2);
    elseif isequal(station_name,'SVALBARD')
        R=diag([deg2rad(125*1e-3) deg2rad(125*1e-3) 0.01].^2);
    else
        error('Incorrect station name')
    end

    % Earth gravitational parameter
    mu=cspice_bodvrd('Earth','GM',1);
    
    % propagation over reference time grid
    if j2==1
        [~,xx_kep]=ode113(@(t,x) KEPJ2(t,x,mu),tt,x0,opt);
    else
        [~,xx_kep]=ode113(@(t,x) KEP(t,x,mu),tt,x0,opt);
    end

    % measurements associated to propagated state
    [exp_meas_kep(1,:),exp_meas_kep(2,:),exp_meas_kep(3,:)]=antenna_pointing(station_name,tt(2:end),xx_kep(2:end,:));

    % residuals matrix and array
    RES=sqrtm(R)\(exp_meas_kep-meas);
    res=RES(:);
end