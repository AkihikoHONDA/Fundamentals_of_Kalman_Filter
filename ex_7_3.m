%% 問題設定
% 物理パラメータ値の設定
rho = 1.23; g = 9.81; eta = 6e3;
M = 3e4; a = 3e4;

% 観測条件
T = 0.5;             % オイラー法離散化周期
EndTime = 30;        % シミュレーション終了時刻
time = 0:T:EndTime;
N = EndTime/T+1;     % データ数
n = 3;               % 状態変数の数
R = 4e3;             % 観測雑音の分散
Q = 0;               % システム雑音はないものとする
B = [0;0;0];

% システム
f = @(x) [x(1)+T*x(2);
    x(2)+T*(0.5*rho*exp(-x(1)/eta)*x(2)^2*x(3)-g);
    x(3)];
h = @(x) sqrt(M^2+(x(1)-a).^2);

% f のヤコビアン(EKFで必要)
A = @(x) [1 T 0;
    -0.5*T*rho/eta*exp(-x(1)/eta)*x(2)^2*x(3) ...
    1+T*rho*exp(-x(1)/eta)*x(2)+x(3) ...
    0.5*T*rho*exp(-x(1)/eta)*x(2)^2;
    0 0 1];

% h のヤコビアン(EKF)で必要
C = @(x) [(x(1)-a)/sqrt(M^2+(x(1)-a)^2); 0; 0];

%% 観測データの生成
w = randn(N,1)*sqrtm(R);  % 観測雑音

% 記憶領域の確保
x = zeros(N,n); y=zeros(n,1);

% 初期条件
x(1,:) = [90000; -6000; 0.03]';
y(1) = h(x(1,:))+w(1);

% 時間更新
for k=2:N
    x(k,:) = f(x(k-1,:));
    y(k) = h(x(k,:))+w(k);
end

%% EKFアルゴリズム
% 推定値記憶領域の確保
xhat_ekf = zeros(N,3);

% 初期値の設定
xhat_ekf (1,:) = x(1,:);             % 状態推定値
P_ekf = [9e3 0 0; 0 4e5 0; 0 0 0.4]; % 誤差共分散

% 推定値を反復更新
for k=2:N
    [xhat_ekf(k,:), P_ekf] = ...
        ekf(f,h,A,B,C,Q,R,y(k),xhat_ekf(k-1,:),P_ekf);
end

%% UKFアルゴリズム
% 推定値記憶領域の確保
xhat_ukf = zeros(N,3);

% 初期値
xhat_ukf(1,:) = x(1,:);              % 状態推定値
P_ukf = [9e3 0 0; 0 4e5 0; 0 0 0.4]; % 誤差共分散

% 推定値を反復更新
for k=2:N
    [xhat_ukf(k,:), P_ukf] = ...
        ukf(f,h,B,Q,R,y(k),xhat_ukf(k-1,:),P_ukf);
end

%% 結果の表示
% 状態の真値の図示
figure(1), clf
for p=1:3
    subplot(3,1,p)
    plot(time,x(:,p));
    xlabel('Time [s]'), ylabel(sprintf('x%d', p))
end

% 推定値の図示
figure(2), clf
for p=1:3
    subplot(3,1,p)
    plot(time, x(:,p), 'k', ...
        time, xhat_ekf(:,p), 'b-.', ...
        time, xhat_ukf(:,p), 'r:');
    xlabel('Time [s]'), ylabel(sprintf('x%d', p))
legend('true', 'ekf', 'ukf')
end

% RMSEの表示
RMSE = @(x) sqrt(mean(x.^2));
fprintf('%10s %10s %10s\n', 'variable', 'RMSE(ekf)', 'RMSE(ukf)');
for p=1:3
    vname = sprintf('x%d', p);
    fprintf('%10s % 10.5f % 10.5f\n', ...
        vname, RMSE(xhat_ekf(:,p)-x(:,p)), ...
        RMSE(xhat_ukf(:,p)-x(:,p)));
end

