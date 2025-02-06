%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 开环系统绘图代码 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% Code for 2D MJS H∞ Asynchronous Controller 双模异步控制器状态轨迹绘制 %%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; % 清空命令行 工作区
tic  % 计时器
system_parameters_tnse % 导入系统参数
%% 绘图参数
T = 80; dt = 1; Nt = T/dt; % 参数：水平方向
Lx = 80; dx = 1; Nx = Lx/dx; % 参数：垂直方向
xh = zeros(Nt, Nx); % xh系统状态信号存储数组
xv = zeros(Nt, Nx); % xv系统状态信号存储数组
w = zeros(Nt, Nx);  % 噪声信号存储数组
%% 水平方向边界条件 initail conditions
for j = 1 : Nx
    if j>=1 && j<=10
        xh(1,j) = 0;
    else
        xh(1,j) = 0;
    end
end
%% 垂直方向边界条件 initail conditions
for i = 1 : Nt 
    if i>=1 && i<=10
        xv(i,1) = 0;
    else
        xv(i,1) = 0;
    end
end
%% 噪声信号
for i = 1 : Nt
    for j = 1 : Nx
        w(i,j) = cos(0.314*(i + j))*exp(-0.15*(i + j));
    end
end
%% 系统与异步滤波器模态切换代码
SystemSeq = zeros(Nt+1,Nx+1);   % 系统状态序列存储数组
ControllerSeq = zeros(Nt+1,Nx+1);   % 异步控制器状态序列存储数组
% flag = round(rand+1);   % 随机生成一个系统模态，1或者2
flag = 1;
for i = 1:Nt+1
    for j = 1:Nx+1
        if flag == 1
            a = rand;  b = rand;
                if a < fiao(1,1)        
                    flagSystem = 1;       % System jump to mode 1
                        if b < fiac(1,1)
                            flagController = 1;
                        else
                            flagController = 2;
                        end
                else      
                    flagSystem = 2;       % System jump to mode 2
                        if b < fiac(2,1)
                            flagController = 1;
                        else
                            flagController = 2;
                        end
                end
        else
            a = rand;  b = rand;
                if a < fiao(2,1)        
                    flagSystem = 1;       % System jump to mode 1
                        if b < fiac(1,1)
                            flagController = 1;
                        else
                            flagController = 2;
                        end
                else      
                    flagSystem = 2;       % System jump to mode 2
                        if b < fiac(2,1)
                            flagController = 1;
                        else
                            flagController = 2;
                        end
                end
        end
        SystemSeq(i,j) = flagSystem;  ControllerSeq(i,j) = flagController;
    end
end
%% 绘图
for i = 1:Nt
    for j = 1:Nx
        if SystemSeq(i,j) == 1 % 系统模式切换
            As = A1; Bs = B1; Es = E1; Cs = C1; Ds = D1; Fs = F1;
        else
            As = A2; Bs = B2; Es = E2; Cs = C2; Ds = D2; Fs = F2;
        end
        % 系统的动态方程
        xh(i+1,j) = As(1,1)*xh(i,j)+As(1,2)*xv(i,j)+Es(1,1)*w(i,j);
        xv(i,j+1) = As(2,1)*xh(i,j)+As(2,2)*xv(i,j)+Es(2,1)*w(i,j);
    end
end
%%
load('xh_state_seqs.mat');  load('xv_state_seqs.mat');  % 导入数据
%% 垂直方向系统状态图
figure('name', '垂直方向系统开环演化图')
[x,y] = meshgrid(1:dt:T, 1:dx:Lx);
z = xv(1:Nt,1:Nx);  mesh(x, y, z);  colormap(winter);
set(gca,'FontSize',15, 'Linewidth', 1);
xlabel('$t$','Interpreter','latex','Fontsize',20);
ylabel('$s$','Interpreter','latex','Fontsize',20);
zlabel('$x^{v}_{s,t}$', 'Interpreter', 'latex', 'Fontsize', 18);
axis([0 T 0 Lx]);  view(-37.8,34.9);  grid on;  grid minor;  box on;
%% 水平方向系统状态图
figure('name', '水平方向系统开环演化图')
[x,y] = meshgrid(1:dt:T, 1:dx:Lx);
z = xh(1:Nt,1:Nx);  mesh(x, y, z);  colormap(winter);
set(gca,'FontSize',15, 'Linewidth', 1);
xlabel('$t$','Interpreter','latex','Fontsize',20);
ylabel('$s$','Interpreter','latex','Fontsize',20);
zlabel('$x^{h}_{s,t}$', 'Interpreter', 'latex', 'Fontsize', 18);
axis([0 T 0 Lx]);  view(-37.8,34.9);  grid on;  grid minor;  box on;
disp(['运行时间: ',num2str(toc)]); 