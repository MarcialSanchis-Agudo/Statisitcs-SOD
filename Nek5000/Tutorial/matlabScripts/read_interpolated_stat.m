%%

clc
close all
clear
format long
double precision;

addpath matFunctions/

%% Setup info

PATH_MESH='../run_interp3d/ZSTAT/';
PATH_INTa='../run_interp3d/output_a_full_30/';
PATH_INTb='../run_interp3d/output_b_full_30/';

Re=10000;
nu = 1/Re;
rho = 1;

%% Load Data

disp('load data')

% Read interpolated data (a files)

stat_a = load_stat_a(PATH_INTa);

% Read interpolated data (b files)

stat_b = load_stat_b(PATH_INTb);

disp('done')

run matFunctions/check_stats.m

%% complete stats


stat_a.uu = stat_b.uu;
stat_a.uw = stat_b.uw;
stat_a.ppp = stat_b.ppp;
stat_a.www = stat_b.www;
stat_a.vvw = stat_b.vvw;
stat_a.Pxx = stat_b.Pxx;
stat_a.Pxz = stat_b.Pxz;
stat_a.Dzz = stat_b.Dzz;
stat_a.Txx = stat_b.Txx;
stat_a.Txz = stat_b.Txz;
stat_a.VDzz = stat_b.VDzz;
stat_a.Pixx = stat_b.Pixx;
stat_a.Pixz = stat_b.Pixz;
stat_a.Czz = stat_b.Czz;
stat_a.Pk = stat_b.Pk;
stat_a.Pik = stat_b.Pik;
stat_a.PTyy = stat_b.PTyy;
stat_a.PTyz = stat_b.PTyz;
stat_a.PSxy = stat_b.PSxy;
stat_a.dUdy = stat_b.dUdy;
stat_a.dVdz = stat_b.dVdz;

stat_a.dPdx = stat_b.dPdx;
stat_a.dPdy = stat_b.dPdy;
stat_a.dPdz = stat_b.dPdz;

%%

save('tut.mat','stat_a','stat_b')
% save('tutHR.mat','stat_a','stat_b')


%% figures

step=0.01;

xx=-10:step:6;
yy=0:step:2.25;

Nx = length(xx);
Ny = length(yy);

[XX,YY]=meshgrid(xx,yy);


U= reshape(stat_a.U,[Ny,Nx]);
V= reshape(stat_a.V,[Ny, Nx]);
W= reshape(stat_a.W,[Ny, Nx]);

dUdy= reshape(stat_a.dUdy,[Ny, Nx]);
tau_wall = (nu*rho)*dUdy(1,:);
u_tau = sqrt(tau_wall/rho);

U_ref =1;
for i=1:Nx
    delta_star(i) = trapz(yy, 1 - (U(:,i)/U_ref));
    a = U(:,i) >= 0.99*U_ref;
    delta_99(i) = min(yy(a));
    theta(i) = trapz(yy, (U(:,i)/U_ref).*(1 - U(:,i)/U_ref));
end

Re_tau = u_tau.*delta_99/nu;
Re_theta = U_ref.*theta/nu;
figure();
plot(xx,delta_star);
grid on;
title('$\delta^*(x)$','FontSize',16,'Interpreter','latex')
xlabel('$x$','FontSize',16,'Interpreter','latex')
ylabel('$\delta^*$','FontSize',16,'Interpreter','latex')

figure();
plot(xx,delta_99);
grid on;
title('$\delta_{99}(x)$','FontSize',16,'Interpreter','latex')
xlabel('$x$','FontSize',16,'Interpreter','latex')
ylabel('$\delta_{99}$','FontSize',16,'Interpreter','latex')

figure()
plot(xx, theta);
grid on;
title('$\theta(x)$','FontSize',16,'Interpreter','latex')
xlabel('$x$','FontSize',16,'Interpreter','latex')
ylabel('$\theta$','FontSize',16,'Interpreter','latex')

figure()
plot(xx, u_tau);
grid on;
title('$u_{\tau}(x)$','FontSize',16,'Interpreter','latex')
xlabel('$x$','FontSize',16,'Interpreter','latex')
ylabel('$u_{\tau}$','FontSize',16,'Interpreter','latex')

figure();
plot(xx,Re_tau);
grid on;
title('$Re_{\tau}(x)$','FontSize',16,'Interpreter','latex')
xlabel('$x$','FontSize',16,'Interpreter','latex')
ylabel('$Re_{\tau}$','FontSize',16,'Interpreter','latex')

figure();
plot(xx,Re_theta);
grid on;
title('$Re_{\theta}(x)$','FontSize',16,'Interpreter','latex')
xlabel('$x$','FontSize',16,'Interpreter','latex')
ylabel('$Re_{\theta}$','FontSize',16,'Interpreter','latex')

% uu= reshape(stat_a.uu,[Ny, Nx]);
% vv= reshape(stat_a.vv,[Ny, Nx]);
% ww= reshape(stat_a.ww,[Ny, Nx]);
% uv= reshape(stat_a.uv,[Ny, Nx]);
% uw= reshape(stat_a.uw,[Ny, Nx]);
% vw= reshape(stat_a.vw,[Ny, Nx]);
% 
% figure();
% contourf(XX, YY, uu);

%%

figure
h=pcolor(XX,YY,U); set(h, 'EdgeColor', 'none'); colorbar;
axis equal tight
title('$U$','FontSize',16,'Interpreter','latex')
xlabel('$x$','FontSize',16,'Interpreter','latex')
ylabel('$y$','FontSize',16,'Interpreter','latex')
% 
% %%
% 
% figure('rend','painters','pos',[10 10 1500 600])
% 
% subplot(1,3,1)
% h=pcolor(XX,YY,U); set(h, 'EdgeColor', 'none'); colorbar;
% axis equal tight
% title('$U$','FontSize',16,'Interpreter','latex')
% xlabel('$x$','FontSize',16,'Interpreter','latex')
% ylabel('$y$','FontSize',16,'Interpreter','latex')
% 
% subplot(1,3,2)
% h=pcolor(XX,YY,V); set(h, 'EdgeColor', 'none'); colorbar;
% axis equal tight
% title('$V$','FontSize',16,'Interpreter','latex')
% xlabel('$x$','FontSize',16,'Interpreter','latex')
% ylabel('$y$','FontSize',16,'Interpreter','latex')
% 
% subplot(1,3,3)
% h=pcolor(XX,YY,W); set(h, 'EdgeColor', 'none'); colorbar;
% axis equal tight
% title('$W$','FontSize',16,'Interpreter','latex')
% xlabel('$x$','FontSize',16,'Interpreter','latex')
% ylabel('$y$','FontSize',16,'Interpreter','latex')
% 
% saveas(gcf,'U.png')
% 
% % 
% 
% uu= reshape(stat_a.uu,[Ny, Nx]);
% vv= reshape(stat_a.vv,[Ny, Nx]);
% ww= reshape(stat_a.ww,[Ny, Nx]);
% uv= reshape(stat_a.uv,[Ny, Nx]);
% uw= reshape(stat_a.uw,[Ny, Nx]);
% vw= reshape(stat_a.vw,[Ny, Nx]);
% 
% figure('rend','painters','pos',[10 10 1600 800])
% 
% subplot(2,3,1)
% h=pcolor(XX,YY,uu); set(h, 'EdgeColor', 'none'); colorbar;
% axis equal tight
% title('$\overline{uu}$','FontSize',16,'Interpreter','latex')
% xlabel('$x$','FontSize',16,'Interpreter','latex')
% ylabel('$y$','FontSize',16,'Interpreter','latex')
% 
% subplot(2,3,2)
% h=pcolor(XX,YY,vv); set(h, 'EdgeColor', 'none'); colorbar;
% axis equal tight
% title('$\overline{vv}$','FontSize',16,'Interpreter','latex')
% xlabel('$x$','FontSize',16,'Interpreter','latex')
% ylabel('$y$','FontSize',16,'Interpreter','latex')
% 
% subplot(2,3,3)
% h=pcolor(XX,YY,ww); set(h, 'EdgeColor', 'none'); colorbar;
% axis equal tight
% title('$\overline{ww}$','FontSize',16,'Interpreter','latex')
% xlabel('$x$','FontSize',16,'Interpreter','latex')
% ylabel('$y$','FontSize',16,'Interpreter','latex')
% 
% subplot(2,3,4)
% h=pcolor(XX,YY,uv); set(h, 'EdgeColor', 'none'); colorbar;
% axis equal tight
% title('$\overline{uv}$','FontSize',16,'Interpreter','latex')
% xlabel('$x$','FontSize',16,'Interpreter','latex')
% ylabel('$y$','FontSize',16,'Interpreter','latex')
% 
% subplot(2,3,5)
% h=pcolor(XX,YY,uw); set(h, 'EdgeColor', 'none'); colorbar;
% axis equal tight
% title('$\overline{uw}$','FontSize',16,'Interpreter','latex')
% xlabel('$x$','FontSize',16,'Interpreter','latex')
% ylabel('$y$','FontSize',16,'Interpreter','latex')
% 
% subplot(2,3,6)
% h=pcolor(XX,YY,vw); set(h, 'EdgeColor', 'none'); colorbar;
% axis equal tight
% title('$\overline{vw}$','FontSize',16,'Interpreter','latex')
% xlabel('$x$','FontSize',16,'Interpreter','latex')
% ylabel('$y$','FontSize',16,'Interpreter','latex')
% 
% saveas(gcf,'uu.png')

%% Duct convergenace

dPdx = reshape(stat_a.dPdx,[Ny, Nx]);

dUdy = reshape(stat_a.dUdy, [Ny, Nx]);
d2Udy2 = diff(dUdy,[], 1);

uv = reshape(stat_a.uv, [Ny,Nx]);
duvdy = diff(uv, [], 1);

V = reshape(stat_a.V, [Ny, Nx]);

dUdz = reshape(stat_a.dUdz, [Ny, Nx]);

index_centerplane_y = round(Ny/2);
index_centerplane_x = round(Nx/2);

eps_duct = -dPdx(index_centerplane_y, index_centerplane_x)...
    + (1/Re)*d2Udy2(index_centerplane_y, index_centerplane_x)...
    - duvdy(index_centerplane_y, index_centerplane_x)...
    - (...
    V(index_centerplane_y, index_centerplane_x) ...
    * dUdy(index_centerplane_y, index_centerplane_x)...
    )...
    + (1/Re)*0;