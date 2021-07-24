function [xhat_new, P_new, G] = ukf(f,h,B,Q,R,y,xhat,P)
% UKFの更新式
%
% [xhat_new, P_new, G] = ukf(f,h,B,Q,R,y,xhat,P)
% UKFの推定値更新を行う
% 引数：
%     f,h,B: 対象システム
%               x(k+1) = f(x(k)) + Bv(k)
%                 y(k) = h(x(k)) + w(k)
%            を記述する関数への関数ハンドル f, h および行列 B
%     注意：対象システムが既知の制御入力 u を持つ関数 fu(x(k),u(k))
%           で記述される場合
%              f = @(x) fu(x,u(k))
%           を与えればよい．
%     Q,R: 雑音v,wの共分散行列．v,w は正規性白色雑音で
%                E[v(k)] = E[w(k)] = 0
%           E[v(k)'v(k)] = Q, E[w(k)'w(k)] = R
%          であることを想定
%     y: 状態更新後時点での観測出力 y(k)
%     xhat,P: 更新前の状態推定値 xhat(k-1)・誤差共分散行列 P(k-1)
%
% 戻り値
%     xhat_new: 更新後の状態推定値 xhat(k)
%     P_new:    更新後の誤差共分散行列 P(k)
%     G:        カルマンゲイン G(k)
%
% 参考：
%     線形カルマンフィルタ: KF
%     拡張カルマンフィルタ: EKF

% 列ベクトルに整形
xhat=xhat(:); y=y(:);

% 事前推定値
[xhatm,Pm] = ut(f,xhat,P);        % U変換による遷移後状態の近似
Pm         = Pm+ B*Q*B';          % システム雑音を考慮
[yhatm,Pyy,Pxy] = ut(h,xhatm,Pm); % U変換による出力値の近似

% カルマンゲイン行列
G = Pxy/(Pyy+R);

% 事後推定値
xhat_new = xhatm + G*(y-yhatm);   % 状態
P_new    = Pm - G*Pxy';           % 誤差共分散
end