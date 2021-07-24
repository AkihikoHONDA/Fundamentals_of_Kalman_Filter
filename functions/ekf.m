function [xhat_new, P_new, G] = ekf(f,h,A,B,C,Q,R,y,xhat,P)
% EKF 拡張カルマンフィルタの更新式
%
% [xhat_new, P_new, G] = ekf(f,h,A,B,C,Q,R,t,xhat,P)
% 拡張カルマンフィルタの推定値更新を行う
% 引数：
%     f,h,B: 対象システム
%               x(k+1) = f(x(k)) + Bv(k)
%                 y(k) = h(x(k)) + w(k)
%            を記述する関数への関数ハンドル f, h および行列 B
%     注意：対象システムが既知の制御入力 u を持つ関数 fu(x(k),u(k))
%           で記述される場合
%              f = @(x) fu(x,u(k))
%           を与えればよい．
%     A,C: f,hのヤコビアンを計算する関数への関数ハンドル
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
%     Unscented カルマンフィルタ: UKF

% 列ベクトルに整形
xhat=xhat(:); y=y(:);

% 事前推定値
xhatm = f(xhat);                             % 状態
Pm = A(xhat)*P*A(xhat)' + B*Q*B';            % 誤差共分散

% カルマンゲイン
G = Pm*C(xhatm)/(C(xhatm)'*Pm*C(xhatm)+R);

% 事後推定値
xhat_new = xhatm+G*(y-h(xhatm));             % 状態
P_new = (eye(size(A(xhat)))-G*C(xhatm)')*Pm; % 誤差共分散
end

