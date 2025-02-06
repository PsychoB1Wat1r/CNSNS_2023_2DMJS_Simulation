
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%% 闭环系统绘图代码 + 2023-9-4 %%%%%%%%%%%%%%%%%%%%%%
%% %%%%%% 2D MJS H∞ Asynchronous Controller 双模异步控制器状态轨迹绘制 %%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; % 清空命令行 工作区
tic  % 计时器
lmi_2detc_cnsns % 调用 LMIs 数值求解代码
%% 绘图坐标 绘图线条 坐标轴字体 大小
fontSizeXY = 15;  lineWidth = 1;  fontSizeAxis = 20; markerSize = 9;
%% 触发次数初始化/数据减少率计算
Trigger_times_h = 0;  Trigger_times_v = 0; % 触发次数计数器初始化
Trigger_times_h_sum = 0;  Trigger_times_v_sum = 0; % 重复模拟 触发次数计数器初始化
Num_of_Repeat = 1; % 重复模拟次数 Monte Carlo
%% 绘图参数
rand('state',10);
T = 30; dt = 1; Nt = T/dt;  % 参数：水平方向
Lx = 30; dx = 1; Nx = Lx/dx;    % 参数：垂直方向
xh = zeros(Nt, Nx); % xh 系统状态信号存储数组
xv = zeros(Nt, Nx); % xv 系统状态信号存储数组
w = zeros(Nt, Nx);  % w 噪声信号存储数组
zh = zeros(Nt, Nx); % zh 控制输出信号存储数组
zv = zeros(Nt, Nx); % zv 控制输出信号存储数组
eventh_array = zeros(Nt, Nx);    % xh 触发瞬时存储数组
eventv_array = zeros(Nt, Nx);    % xv 触发瞬时存储数组
eventh_time_array = zeros(Nt, Nx);    % 保存水平事件时间以计算采样时间
eventv_time_array = zeros(Nt, Nx);    % 保存水平事件时间以计算采样时间
zhsum = zeros(Nt, Nx);  % zh H无穷性能指标演化存储数组
zvsum = zeros(Nt, Nx);  % zv H无穷性能指标演化存储数组
PackerlossSeqh = -ones(Nt, Nx); % 水平方向信号丢包过程
PackerlossSeqv = -ones(Nt, Nx); % 垂直方向信号丢包过程
% 
u = 0;  u_old = 0;
u_array = zeros(Nt, Nx); % u 控制信号存储数组
uh = 0;  uh_old = 0;
uh_array = zeros(Nt, Nx); % uh 控制信号存储数组
uv = 0;  uv_old = 0;
uv_array = zeros(Nt, Nx); % uv 控制信号存储数组
xh_new = zeros(1, 1);  xv_new = zeros(1, 1); % 触发控制初始时刻系统状态值选取为零 current
xh_old = zeros(1, 1);  xv_old = zeros(1, 1); % 触发控制初始时刻系统状态值选取为零 last
xhiter = xh(1,1);  xviter = xv(1,1); % 迭代过程的中间变量（迭代器）
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
%% 系统与异步控制器模态切换代码
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
        SystemSeq(i,j) = flagSystem; 
        ControllerSeq(i,j) = flagController;
    end
end
%% 绘图/减少率计算
for r = 1 : Num_of_Repeat % 重复模拟
    Trigger_times_h = 0;  Trigger_times_v = 0; % 触发次数计数器
    for i = 1 : Nt
        last_eventh = 0;  last_eventv = 0; % 水平与垂直方向上 上一次事件出现的时刻
        for j = 1 : Nx
            if  ControllerSeq(i,j) == 1 % 异步控制器模态的切换
                K = K1; % 选择异步控制器一
            else
                K = K2; % 选择异步控制器二
            end
            if SystemSeq(i,j) == 1  % 系统模态切换
                As = A1;  Bs = B1;  Es = E1;  Cs = C1;  Ds = D1;  Fs = F1;  Omega = Omega1;
            else
                As = A2;  Bs = B2;  Es = E2;  Cs = C2;  Ds = D2;  Fs = F2;  Omega = Omega2;
            end
            xh_error = xh(i,j) - xh_old;  xv_error = xv(i,j) - xv_old; % 初始化状态误差
           %% 事件触发生成器
            % 水平方向事件触发生成器
            if ( Omega(1,1)*abs(norm(xh_error)) >= etc(1,1)*Omega(1,1)*(norm(xh(i,j))) )   % 水平触发条件成立
                eventh_array(i,j) = 1;  % 用于绘制水平触发时刻分布图
                eventh_time_array(i,j) = j - last_eventh;  % 用于绘制触发时刻的分布图
                last_eventh = j;
                Trigger_times_h = Trigger_times_h + 1; % 水平触发次数加1
            else  
                eventh_array(i,j) = 0; % 用于绘制水平触发时刻分布图
                eventh_time_array(i,j) = 0;
            end
            % 垂直方向事件触发生成器
            if ( Omega(2,2)*abs(norm(xv_error)) >= etc(2,2)*Omega(2,2)*(norm(xv(i,j))) )   % 垂直触发条件成立
                eventv_array(i,j) = 1; % 用于绘制垂直触发时刻分布图
                eventv_time_array(i,j) = j - last_eventv;  % 用于绘制触发时刻的分布图
                last_eventv = j;
                Trigger_times_v = Trigger_times_v + 1; % 垂直触发次数加1
            else  
                eventv_array(i,j) = 0; % 用于绘制垂直触发时刻 分布图
                eventv_time_array(i,j) = 0;
            end
            %% 测量信号更新
            if eventh_array(i,j) == 1
                probh = rand; % 生成一个随机概率用于判断是否发生丢包行为
                if probh <= k
                    flagPacketh = 1; % 不丢包
                    xh_new = xh(i,j); % xh 信号更新传输
                    xh_old = xh_new; % 最后一次触发时刻的系统状态更新
                else
                    flagPacketh = 0; % 丢包
                    xh_new = 0; 
                    xh_old = xh_new; % 最后一次触发时刻的系统状态更新
                end
                PackerlossSeqh(i,j) = flagPacketh; % 存储丢包序列 用于仿真模拟
            else
                xh_new = xh_old;
            end
            if eventv_array(i,j) == 1
                probv = rand; % 生成一个随机概率用于判断是否发生丢包行为
                if probv <= k
                    flagPacketv = 1; % 不丢包
                    xv_new = xv(i,j); % xv 信号更新传输
                    xv_old = xv_new; % 最后一次触发时刻的系统状态更新
                else
                    flagPacketv = 0; % 丢包
                    xh_new = 0; 
                    xv_old = xv_new; % 最后一次触发时刻的系统状态更新
                end
                PackerlossSeqv(i,j) = flagPacketv; % 存储丢包序列 用于仿真模拟
            else
                xv_new = xv_old;
            end
            u = K*[xh_new  xv_new]'; % 控制信号更新
            u_array(i,j) = u; % 存储控制信号 u 用于绘图
            %% 系统响应
            xh(i+1,j) = As(1,1)*xh(i,j) + As(1,2)*xv(i,j) + Bs(1,:)*u + Es(1,1)*w(i,j);
            xv(i,j+1) = As(2,1)*xh(i,j) + As(2,2)*xv(i,j) + Bs(2,:)*u + Es(2,1)*w(i,j);
            zh(i+1,j) = Cs(1,1)*xh(i,j) + Cs(1,2)*xv(i,j) + Ds(1,:)*u + Fs(1,1)*w(i,j);
            zv(i,j+1) = Cs(2,1)*xh(i,j) + Cs(2,2)*xv(i,j) + Ds(2,:)*u + Fs(2,1)*w(i,j);
            
            zhsum(i,j) = sqrt(norm(zh(1:i,1:j), 2))/sqrt(norm(w(1:i,1:j), 2)); % zh 用于绘制性能指标函数图像（范数）
            zvsum(i,j) = sqrt(norm(zv(1:i,1:j), 2))/sqrt(norm(w(1:i,1:j), 2)); % zv 用于绘制性能指标函数图像（范数）
        end 
    end
    Trigger_times_h_sum = Trigger_times_h_sum + Trigger_times_h; % 水平触发次数 均值求解 暂存变量
    Trigger_times_v_sum = Trigger_times_v_sum + Trigger_times_v; % 垂直触发次数 均值求解 暂存变量
end
%% 垂直方向系统状态图
figure('name', '垂直方向系统闭环演化图')
[x,y] = meshgrid(2:dt:T, 2:dx:Lx);
z = xv(2:dt:Nt, 2:dx:Nx);
mesh(x, y, z, 'LineWidth', lineWidth);
set(gca,'FontSize', fontSizeXY, 'Linewidth', lineWidth);
xlabel('$t$', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis);
ylabel('$s$', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis);
zlabel('$x^{v}_{s,t}$', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis-2);
axis([0 T 0 Lx]);  view(-37.8,34.9);  colormap(winter);  grid on;  grid minor;  box on;
%% 水平方向系统状态图
figure('name', '水平方向系统闭环演化图')
[x,y] = meshgrid(2:dt:T, 2:dx:Lx);
z = xh(2:dt:Nt, 2:dx:Nx);
mesh(x, y, z, 'LineWidth', lineWidth);
set(gca, 'FontSize', fontSizeXY, 'Linewidth', lineWidth);
xlabel('$t$', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis);
ylabel('$s$', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis);
zlabel('$x^{h}_{s,t}$', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis-2);
axis([0 T 0 Lx]);  view(-37.8,34.9);  colormap(winter);  grid on;  grid minor;  box on;
%% 异步控制输入信号
figure('name', '控制输入演化图')
[x,y] = meshgrid(2:dt:T, 2:dx:Lx);
z = u_array(2:dt:Nt, 2:dx:Nx);
mesh(x, y, z, 'LineWidth', lineWidth);
set(gca, 'FontSize', fontSizeXY, 'Linewidth', lineWidth);
xlabel('$t$', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis);
ylabel('$s$', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis);
zlabel('$u_{s,t}$', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis-2);
axis([0 T 0 Lx]);  view(-37.8,34.9);  colormap(winter);  grid on;  grid minor;  box on;
%% zh 性能比较
figure('name', 'zh 性能比较');
[x,y] = meshgrid(1:dt:T, 1:dx:Lx);
z = zhsum(1:dt:T, 1:dt:Lx);
mesh(x, y, z, 'LineWidth', lineWidth);
set(gca, 'FontSize', fontSizeXY, 'Linewidth', lineWidth);
xlabel('$\hat{t}$', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis);
ylabel('$\hat{s}$', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis);
zlabel('$\gamma^{h}_{\hat{s},\hat{t}}$', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis-2);
axis([0 T 0 Lx 0 0.6]);  view(-37.8,34.9);  colormap(hsv(3));  grid on;  grid minor;  box on;
%% zv 性能比较
figure('name', 'zv 性能比较');
[x,y] = meshgrid(1:dt:T, 1:dx:Lx);
z = zvsum(1:dt:T, 1:dt:Lx);
mesh(x, y, z, 'LineWidth', lineWidth);
set(gca, 'FontSize', fontSizeXY, 'Linewidth', lineWidth);
xlabel('$\hat{t}$', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis);
ylabel('$\hat{s}$', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis);
zlabel('$\gamma^{v}_{\hat{s},\hat{t}}$', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis-2);
axis([0 T 0 Lx 0 0.6]);  view(-37.8,34.9);  colormap(hsv(3));  grid on;  grid minor;  box on;
%% 系统模式切换
figure('name', '系统模式切换')
[x,y] = meshgrid(0:dt:T, 0:dx:Lx);
z = SystemSeq(1:dt:T+1, 1:dt:Lx+1);
surf(x,y,z);  view(0,90);  colormap(cool);
set(gca, 'FontSize', fontSizeXY, 'Linewidth', lineWidth);
xlabel('$t$', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis);
ylabel('$s$', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis);
axis([0 T 0 Lx 0.5 2.5]);  grid on;  grid minor;  box on;
%% 控制器模式切换
figure('name', '控制器模式切换')
[x,y] = meshgrid(0:dt:T, 0:dx:Lx);
z = ControllerSeq(1:dt:T+1,1:dt:Lx+1);
surf(x,y,z);  view(0,90);  colormap(cool);
set(gca,'FontSize',fontSizeXY, 'Linewidth', lineWidth);
xlabel('$t$','Interpreter','latex','Fontsize',fontSizeAxis);
ylabel('$s$','Interpreter','latex','Fontsize',fontSizeAxis);
axis([0 T 0 Lx 0.5 2.5]);  grid on;  grid minor;  box on;
%% 事件触发释放时刻 xh
figure('name', '事件触发释放时刻 xh');
h = stem3((0:dt:T-1), (0:dt:T-1), eventh_time_array(1:dt:T, 1:dt:Lx), 'linewidth', lineWidth/2, 'marker', '.', 'color', '#1D267D', 'markersize', 13);
xdata = get(h,'xdata');  ydata = get(h,'ydata');  zdata = get(h,'zdata');   % 通过句柄获取绘图数据的 x y z 坐标值
zdata(zdata == 0) = NaN;    % 将 z = 0 的值重新置为 NaN
set(h,'xdata', xdata,'ydata', ydata, 'zdata', zdata);     % 更新绘图数据 (x 和 y 坐标值)
set(gca,'FontSize', fontSizeXY, 'Linewidth', lineWidth);
xlabel('$t$', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis);
ylabel('$s$', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis);
zlabel('Release interval', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis-4);
view(-37.8,34.9);  grid on;  grid minor;  box on;
%% 事件触发释放时刻 xv #0A4D68
figure('name', '事件触发释放时刻 xv');
h = stem3((0:dt:T-1), (0:dt:T-1), eventv_time_array(1:dt:T, 1:dt:Lx), 'linewidth', lineWidth/2, 'marker', '.', 'color', '#1D267D', 'markersize', 13);
xdata = get(h,'xdata');  ydata = get(h,'ydata');  zdata = get(h,'zdata');   % 通过句柄获取绘图数据的 x y z 坐标值
zdata(zdata == 0) = NaN;    % 将 z = 0 的值重新置为 NaN
set(h,'xdata', xdata,'ydata', ydata, 'zdata', zdata);     % 更新绘图数据 (x 和 y 坐标值)
set(gca,'FontSize', fontSizeXY, 'Linewidth', lineWidth);
xlabel('$t$', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis);
ylabel('$s$', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis);
zlabel('Release interval', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis-4);
view(-37.8,34.9);  grid on;  grid minor;  box on;
%% 丢包情况 3D 水平方向
figure('name', '水平方向 丢包情况')
h = stem3((0:dt:T-1), (0:dt:T-1), PackerlossSeqh(1:dt:T, 1:dt:Lx), 'linestyle', 'none', 'linewidth', lineWidth, 'marker', '.', 'color', '#0A4D68', 'markersize', markerSize);
xdata = get(h,'xdata');  ydata = get(h,'ydata');  zdata = get(h,'zdata');   % 通过句柄获取绘图数据的 x y z 坐标值
zdata(zdata == -1) = NaN;    % 将 z = -1 的值重新置为 NaN
zdata(zdata == 1) = NaN;    % 将 z = 1 的值重新置为 NaN
set(h,'xdata', xdata, 'ydata', ydata, 'zdata', zdata, 'linestyle', 'none', 'linewidth', lineWidth, 'marker', 'x', 'color', 'r', 'markersize', 6);
hold on 
h1 = stem3((0:dt:T-1), (0:dt:T-1), PackerlossSeqh(1:dt:T, 1:dt:Lx), 'linestyle', 'none', 'linewidth', lineWidth, 'marker', '.', 'color', '#0A4D68', 'markersize', markerSize);
xdata = get(h1,'xdata');  ydata = get(h1,'ydata');  zdata = get(h1,'zdata');   % 通过句柄获取绘图数据的 x y z 坐标值
zdata(zdata == -1) = NaN;    % 将 z = -1 的值重新置为 NaN
zdata(zdata == 0) = NaN;    % 将 z = 0 的值重新置为 NaN
set(h1,'xdata', xdata, 'ydata', ydata, 'zdata', zdata, 'linestyle', 'none', 'linewidth', lineWidth, 'marker', '.', 'color', 'b', 'markersize', markerSize);
set(gca,'FontSize', fontSizeXY, 'Linewidth', lineWidth);
xlabel('$t$', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis);
ylabel('$s$', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis);
h2 = legend('Packet loss', 'Packet arrival', 'Interpreter', 'latex', 'Fontsize', 12, 'location', 'northeast');
% set(h2,'Orientation', 'horizon', 'Box', 'on');
view(0,90);  grid on;  grid minor;  box on;
%% 丢包情况 3D 垂直方向
figure('name', '垂直方向 丢包情况')
h = stem3((0:dt:T-1), (0:dt:T-1), PackerlossSeqv(1:dt:T, 1:dt:Lx), 'linestyle', 'none', 'linewidth', lineWidth, 'marker', '.', 'color', '#0A4D68', 'markersize', markerSize);
xdata = get(h,'xdata');  ydata = get(h,'ydata');  zdata = get(h,'zdata');   % 通过句柄获取绘图数据的 x y z 坐标值
zdata(zdata == -1) = NaN;    % 将 z = -1 的值重新置为 NaN
zdata(zdata == 1) = NaN;    % 将 z = 1 的值重新置为 NaN
set(h,'xdata', xdata, 'ydata', ydata, 'zdata', zdata, 'linestyle', 'none', 'linewidth', lineWidth, 'marker', 'x', 'color', 'r', 'markersize', 6);
hold on 
h1 = stem3((0:dt:T-1), (0:dt:T-1), PackerlossSeqv(1:dt:T, 1:dt:Lx), 'linestyle', 'none', 'linewidth', lineWidth, 'marker', '.', 'color', '#0A4D68', 'markersize', markerSize);
xdata = get(h1,'xdata');  ydata = get(h1,'ydata');  zdata = get(h1,'zdata');   % 通过句柄获取绘图数据的 x y z 坐标值
zdata(zdata == -1) = NaN;    % 将 z = -1 的值重新置为 NaN
zdata(zdata == 0) = NaN;    % 将 z = 0 的值重新置为 NaN
set(h1,'xdata', xdata, 'ydata', ydata, 'zdata', zdata, 'linestyle', 'none', 'linewidth', lineWidth, 'marker', '.', 'color', 'b', 'markersize', markerSize);
set(gca,'FontSize', fontSizeXY, 'Linewidth', lineWidth);
xlabel('$t$', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis);
ylabel('$s$', 'Interpreter', 'latex', 'Fontsize', fontSizeAxis);
h2 = legend('Packet loss', 'Packet arrival', 'Interpreter', 'latex', 'Fontsize', 12, 'location', 'northeast');
% set(h2,'Orientation', 'horizon', 'Box', 'on');
view(0,90);  grid on;  grid minor;  box on;
%% 传输信号减少率
disp('****************************************');
disp('***数据包抵达率/H∞噪声抑制水平/触发参数***');
disp('****************************************');
Data_arrival_rate = k,  Gamma = gamma,  Trigger_thresholds = etc
disp('***********************');
disp('***水平/垂直 触发次数***');
disp('***********************');
Trigger_times_h,  Trigger_times_v
disp('*************************');
disp('***水平/垂直 传输减少率***');
disp('*************************');
Trigger_rate_h = 1 - Trigger_times_h / (Nx*Nt)
Trigger_rate_v = 1 - Trigger_times_v / (Nx*Nt)
disp('****************************');
disp('***水平/垂直 平均传输减少率***');
disp('****************************');
Mean_Trigger_rate_h = 1 - (Trigger_times_h_sum / Num_of_Repeat) / (Nx*Nt)
Mean_Trigger_rate_v = 1 - (Trigger_times_v_sum / Num_of_Repeat) / (Nx*Nt)
disp('****************************');
disp('****水平/垂直 爱趣无穷范数****');
disp('****************************');
Norm_h = sqrt(norm(zh(1:Nt,1:Nx),2)) / sqrt(norm(w(1:Nt,1:Nx), 2))
Norm_v = sqrt(norm(zv(1:Nt,1:Nx),2)) / sqrt(norm(w(1:Nt,1:Nx), 2))
disp(['运行时间: ', num2str(toc)]); 