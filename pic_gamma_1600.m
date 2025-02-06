%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%% 双模态三维性能指标 gamma 绘制 + 2023-9-4 %%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 步长 0.025 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; % 清空命令行 工作区
tic  % 计时器
load('gammaset.mat'); % 导入数据
%% 三维性能指标绘图参数
T = 1; dt = 0.025; Nt = T/dt; % 坐标轴1
Lx = 1; dx = 0.025; Nx = Lx/dx; % 坐标轴2
%% 性能指标gamma 三维图
H1 = figure;
[x,y] = meshgrid(0:dt:T,0:dx:Lx);
z = gammaset(1:Nt+1,1:Nx+1);
mesh(x, y, z, 'LineWidth', 1, 'marker', '.', 'MarkerFaceColor',[0.5,0.5,0.5], 'markersize', 10);
xlabel('$\pi_{21}$','Interpreter','latex','Fontsize',19);
ylabel('$\pi_{11}$','Interpreter','latex','Fontsize',19);
zlabel('$\gamma_{c}^{*}$','Interpreter','latex','Fontsize',19);
axis([0 T 0 Lx 1.52 1.83]);  grid on;  grid minor;  box on;
colormap(hsv(3));
%% 性能指标gamma 二维图
H2 = figure;
plot((0:dt:T), diag(gammaset(1:Nt+1,1:Nx+1)), 'color', '#083AA9', 'linestyle', '-.', 'linewidth', 2);
hold on
plot((0:dt:T), diag(rot90(gammaset(1:Nt+1,1:Nx+1))), 'color', '#FF2E63', 'linestyle', '-.', 'linewidth', 2, 'marker', 'p', 'MarkerFaceColor',[0.5,0.5,0.5], 'markersize', 4.5);
hold on
plot([0 T],[1.53956 1.53956], 'color', '#222831', 'linewidth', 2);
set(gca,'FontSize',11, 'Linewidth', 1);
xlabel('$\pi_{11}$','Interpreter','latex','Fontsize', 19);
ylabel('$\gamma_{c}^{*}$','Interpreter','latex','Fontsize', 19);
h1 = legend('Mode-independent $(\pi_{11}=\pi_{21})$', 'Asynchronous $(\pi_{11}+\pi_{21}=1)$', 'Synchronous $(\pi_{11}=1$ and $\pi_{21}=0)$', 'Interpreter','latex','Fontsize',13, 'location', 'best');
axis([0 T 1.52 1.82]);  grid on;  grid minor;  box on;
disp(['运行时间: ',num2str(toc)]); 