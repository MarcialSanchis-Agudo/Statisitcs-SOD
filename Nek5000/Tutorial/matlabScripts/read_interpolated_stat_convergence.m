%%

clc
close all
clear
format long
double precision;

addpath matFunctions/

%% Setup info

PATH_MESH='../run_interp3d/ZSTAT/';
for i = 1:85
    PATH_INTa=['../run_interp3d/output_a_full_',num2str(i),'/'];
    PATH_INTb=['../run_interp3d/output_b_full_',num2str(i), '/'];


    %% Load Data

    disp('load data')

    % Read interpolated data (a files)

    stat_a(i) = load_stat_a(PATH_INTa);

    % Read interpolated data (b files)

    stat_b(i) = load_stat_b(PATH_INTb);
end

Re=10000;
rho = 1;
nu = 1/Re;



% disp('done')
%
% stat_a = stat_1_a;
% stat_b = stat_1_b;
% disp('Stat 1')
% run matFunctions/check_stats.m
%
% stat_a = stat_2_a;
% stat_b = stat_2_b;
% disp('Stat 2')
% run matFunctions/check_stats.m


%% complete stats

for i = 1:length(stat_a)
    stat_a(i).uu = stat_b(i).uu;
    stat_a(i).uw = stat_b(i).uw;
    stat_a(i).ppp = stat_b(i).ppp;
    stat_a(i).www = stat_b(i).www;
    stat_a(i).vvw = stat_b(i).vvw;
    stat_a(i).Pxx = stat_b(i).Pxx;
    stat_a(i).Pxz = stat_b(i).Pxz;
    stat_a(i).Dzz = stat_b(i).Dzz;
    stat_a(i).Txx = stat_b(i).Txx;
    stat_a(i).Txz = stat_b(i).Txz;
    stat_a(i).VDzz = stat_b(i).VDzz;
    stat_a(i).Pixx = stat_b(i).Pixx;
    stat_a(i).Pixz = stat_b(i).Pixz;
    stat_a(i).Czz = stat_b(i).Czz;
    stat_a(i).Pk = stat_b(i).Pk;
    stat_a(i).Pik = stat_b(i).Pik;
    stat_a(i).PTyy = stat_b(i).PTyy;
    stat_a(i).PTyz = stat_b(i).PTyz;
    stat_a(i).PSxy = stat_b(i).PSxy;
    stat_a(i).dUdy = stat_b(i).dUdy;
    stat_a(i).dVdz = stat_b(i).dVdz;

    stat_a(i).dPdx = stat_b(i).dPdx;
    stat_a(i).dPdy = stat_b(i).dPdy;
    stat_a(i).dPdz = stat_b(i).dPdz;
end

step=0.01;

xx=-10:step:6;
yy=0:step:2.25;
Nx = length(xx);
Ny = length(yy);
x_start = -7;
index_start = find(xx>=x_start,1);

k = zeros(length(stat_a(1).U(index_start*Ny:end)), length(stat_a));
for i = 1:length(stat_a)
    max_U(i) = max(stat_a(i).U(index_start*Ny:end));
    max_V(i) = max(stat_a(i).V(index_start*Ny:end));
    max_W(i) = max(stat_a(i).W(index_start*Ny:end));
    k(:,i) = k(:,i) + ...
        0.5*(...
        stat_a(i).uu(index_start*Ny:end) + ...
        stat_a(i).vv(index_start*Ny:end) + ...
        stat_a(i).ww(index_start*Ny:end)...
        );

    k_total(i) = sum(k(:,i));
end
t_plot = 32:2:200;
figure()
plot(t_plot, k_total)
title('k')
figure()
semilogy(t_plot(1:end), abs(k_total-k_total(end))/k_total(end))
title('k diff')
figure()
plot(t_plot,max_U)
title('max U')
figure()
plot(t_plot, max_V)
title('max V')
figure()
plot(t_plot, max_W)
title('max W')
%%


% %save('tutHR.mat','stat_a','stat_1_b')
%
%
%% figures


[XX,YY]=meshgrid(xx,yy);
index_plot = 85;
U= reshape(stat_a(index_plot).U,[Ny, Nx]);
V= reshape(stat_a(index_plot).V,[Ny, Nx]);
W= reshape(stat_a(index_plot).W,[Ny, Nx]);

uu = reshape(stat_a(index_plot).uu,[Ny, Nx]);
b = uu;
b(b>=0) = nan;


figure
h=pcolor(XX,YY,U); set(h, 'EdgeColor', 'none'); colorbar;
% axis equal tight
title('$U$','FontSize',16,'Interpreter','latex')
xlabel('$x$','FontSize',16,'Interpreter','latex')
ylabel('$y$','FontSize',16,'Interpreter','latex')

figure
h=pcolor(XX,YY,b); set(h, 'EdgeColor', 'none'); colorbar;
% axis equal tight
title('$uu$','FontSize',16,'Interpreter','latex')
xlabel('$x$','FontSize',16,'Interpreter','latex')
ylabel('$y$','FontSize',16,'Interpreter','latex')

dUdy= reshape(stat_a(index_plot).dUdy,[Ny, Nx]);
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



figure();
plot(U(:,xx==-5),yy);
grid on;
hold on;
plot(U(:,xx==0),yy);
plot(U(:,xx==5),yy);
title('$U$','FontSize',16,'Interpreter','latex')
xlabel('$U$','FontSize',16,'Interpreter','latex')
ylabel('$y$','FontSize',16,'Interpreter','latex')

figure();
plot(uu(:,xx==-5),yy);
grid on;
hold on;
plot(uu(:,xx==0),yy);
plot(uu(:,xx==5),yy);
title('$uu$','FontSize',16,'Interpreter','latex')
xlabel('$uu$','FontSize',16,'Interpreter','latex')
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
% %%
%
% uu= reshape(stat_1_a.uu,[Ny, Nx]);
% vv= reshape(stat_1_a.vv,[Ny, Nx]);
% ww= reshape(stat_1_a.ww,[Ny, Nx]);
% uv= reshape(stat_1_a.uv,[Ny, Nx]);
% uw= reshape(stat_1_a.uw,[Ny, Nx]);
% vw= reshape(stat_1_a.vw,[Ny, Nx]);
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

% eps_duct = zeros(1, length(stat_a));
% for i = 1:length(stat_a)
%
%     tau_wall = (nu*rho)*dUdy(1,:);
%     u_tau = sqrt(tau_wall/rho);
%
%     dPdx = reshape(stat_a(i).dPdx,[Ny, Nx]);
%
%     dUdy = reshape(stat_a(i).dUdy, [Ny, Nx]);
%     d2Udy2 = diff(dUdy,[], 1);
%
%     uv = reshape(stat_a(i).uv, [Ny,Nx]);
%     duvdy = diff(uv, [], 1);
%
%     V = reshape(stat_a(i).V, [Ny, Nx]);
%
%     dUdz = reshape(stat_a(i).dUdz, [Ny, Nx]);
%
%     index_centerplane_y = 1:size(d2Udy2,1);
%     index_centerplane_x = round(3*Nx/4);
%
%     eps_duct(i) = rms((...
%             - dPdx(index_centerplane_y, index_centerplane_x)...
%             + (1/Re)*d2Udy2(index_centerplane_y, index_centerplane_x)...
%             - duvdy(index_centerplane_y, index_centerplane_x)...
%             - (...
%                 V(index_centerplane_y, index_centerplane_x) ...
%                 .* dUdy(index_centerplane_y, index_centerplane_x)...
%             )...
%             + (1/Re)*0 ...
%             ).*nu/(u_tau(index_centerplane_x).^3));
% end
% figure();
% plot(1:30, eps_duct)
% grid on;
% title('$\varepsilon_{duct}$','FontSize',16,'Interpreter','latex')
% xlabel('$T$','FontSize',16,'Interpreter','latex')
% ylabel('$\varepsilon_{duct}$','FontSize',16,'Interpreter','latex')

% save('output_full_convergence.mat', '-v7.3')
%% Convergence velocity
% load('output_full_convergence.mat')

U_prev = reshape(stat_a(1).U, [Ny, Nx]);
V_prev = reshape(stat_a(1).V, [Ny, Nx]);
W_prev = reshape(stat_a(1).W, [Ny, Nx]);
uv_prev = reshape(stat_a(1).uv, [Ny,Nx]);
k_prev = reshape(stat_a(1).uu+stat_a(1).vv+stat_a(1).ww, [Ny, Nx]);
resk_prev = reshape(stat_a(1).Resk, [Ny,Nx]);

eps_U_duct = zeros(Ny, Nx, length(stat_a)-1);
eps_V_duct = zeros(Ny, Nx, length(stat_a)-1);
eps_W_duct = zeros(Ny, Nx, length(stat_a)-1);
eps_uv_duct = zeros(Ny, Nx, length(stat_a)-1);

rms_U_domain = zeros(length(stat_a)-1,1);
rms_V_domain = zeros(length(stat_a)-1,1);
rms_W_domain = zeros(length(stat_a)-1,1);
rms_uv_domain = zeros(length(stat_a)-1,1);
rms_k_domain = zeros(length(stat_a)-1,1);

max_U_domain = zeros(length(stat_a)-1,1);
max_V_domain = zeros(length(stat_a)-1,1);
max_W_domain = zeros(length(stat_a)-1,1);
max_uv_domain = zeros(length(stat_a)-1,1);
max_k_domain = zeros(length(stat_a)-1,1);

sum_k_domain = zeros(length(stat_a)-1,1);

max_resk_domain = zeros(length(stat_a)-1,1);
rms_resk_domain = zeros(length(stat_a)-1,1);
l2_resk_domain = zeros(length(stat_a)-1,1);

x_start = -7;
y_end = 0.5;
for i = 2:length(stat_a)

    dUdy = reshape(stat_a(i).dUdy, [Ny, Nx]);
    tau_wall = (nu*rho)*dUdy(1,:);
    u_tau = sqrt(tau_wall/rho);


    U = reshape(stat_a(i).U, [Ny, Nx]);
    V = reshape(stat_a(i).V, [Ny, Nx]);
    W = reshape(stat_a(i).W, [Ny, Nx]);
    uv = reshape(stat_a(i).uv, [Ny,Nx]);
    k = reshape(stat_a(i).uu+stat_a(i).vv+stat_a(i).ww, [Ny, Nx]);
    resk = reshape(stat_a(i).Resk, [Ny,Nx]);

    dUdz = reshape(stat_a(i).dUdz, [Ny, Nx]);


    eps_U_duct(:,:,i-1) = abs(U - U_prev);
    eps_V_duct(:,:,i-1)  = abs(V - V_prev);
    eps_W_duct(:,:,i-1)  = abs(W - W_prev);
    eps_uv_duct(:,:,i-1)  = abs(uv - uv_prev);

    rms_U_domain(i-1) = rms(U(yy<=y_end, xx>=x_start) - U_prev(yy<=y_end, xx>=x_start) , "all");
    rms_V_domain(i-1) = rms(V(yy<=y_end, xx>=x_start) - V_prev(yy<=y_end, xx>=x_start) , "all");
    rms_W_domain(i-1) = rms(W(yy<=y_end, xx>=x_start) - W_prev(yy<=y_end, xx>=x_start) , "all");
    rms_uv_domain(i-1) = rms(uv(yy<=y_end, xx>=x_start) - uv_prev(yy<=y_end, xx>=x_start) , "all");
    rms_k_domain(i-1) = rms(k(yy<=y_end, xx>=x_start) - k_prev(yy<=y_end, xx>=x_start) , "all");

    max_U_domain(i-1)= max(abs(U(yy<=y_end, xx>=x_start) - U_prev(yy<=y_end, xx>=x_start)),[], "all");
    max_V_domain(i-1)= max(abs(V(yy<=y_end, xx>=x_start) - V_prev(yy<=y_end, xx>=x_start)),[], "all");
    max_W_domain(i-1)= max(abs(W(yy<=y_end, xx>=x_start) - W_prev(yy<=y_end, xx>=x_start)),[], "all");
    max_uv_domain(i-1)= max(abs(uv(yy<=y_end, xx>=x_start) - uv_prev(yy<=y_end, xx>=x_start)),[], "all");
    max_k_domain(i-1)= max(abs(k(yy<=y_end, xx>=x_start) - k_prev(yy<=y_end, xx>=x_start)),[], "all");


    sum_k_domain(i-1) = abs(sum(k(yy<=y_end, xx>=x_start),"all") - sum(k_prev(yy<=y_end, xx>=x_start),"all"));

    max_resk_domain(i-1) = max(resk(yy<=y_end, xx>=x_start),[], "all");
    rms_resk_domain(i-1) = rms(resk(yy<=y_end, xx>=x_start), "all");
    l2_resk_domain(i-1) = norm(resk(yy<=y_end, xx>=x_start));



    U_prev = U;
    V_prev = V;
    W_prev = W;
    uv_prev = uv;
    k_prev = k;
    resk_prev = resk;

end

x_plot_index = round(Nx/4);
y_plot_index = round(Ny/10);
t_plot = 34:2:200;

figure();
plot(t_plot, max_resk_domain);
grid on;
title('Convergance $Resk$','FontSize',16,'Interpreter','latex')
xlabel('$T$','FontSize',16,'Interpreter','latex')
ylabel('$max|Resk|$','FontSize',16,'Interpreter','latex')

figure();
plot(t_plot, rms_resk_domain);
grid on;
title('Convergance $Resk$','FontSize',16,'Interpreter','latex')
xlabel('$T$','FontSize',16,'Interpreter','latex')
ylabel('$rms(Resk)$','FontSize',16,'Interpreter','latex')

