function [ym, Pyy, Pxy] = ut(f, xm, Pxx)
% UT(U変換)
% [ym, Pyy, Pxy] = ut(f, xm, Pxx)
%  確率変数xに関して
%    xm   : E[x]
%    Pxx  : E[(x-xm)(x-xm)']
%  が与えられているとき，
%  非線形写像y=f(x)で与えられる確率変数yについて
%    ym   : E[y]
%    Pyy  : E[(y-ym)(y-ym)']
%    Pxy  : E[(x-xm)(y-ym)']
%  をU変換に基づいて計算する．
%  fは関数ハンドルで与えられるものとする．

%% 準備
% 列ベクトルに整形
xm = xm(:);

% mapcols(f,x) : xの各列をfで写像する関数
mapcols = @(f, x) ...
    cell2mat( ...
    cellfun(f, ...
    mat2cell(x, size(x,1), ones(1,size(x,2))) ...
        , 'UniformOutput', false));

% 定数
n     = length(xm);        % 次数
kappa = 3-n;               % スケーリングパラメータ
w0    = kappa/(n + kappa); % 重み
wi    = 1/(2*(n + kappa));
W     = diag([w0; wi*ones(2*n, 1)]);

%% U変換
% シグマポイントの生成
L = chol(Pxx);
X = [xm';
     ones(n,1)*xm'+sqrt(n+kappa)*L;
     ones(n,1)*xm'-sqrt(n+kappa)*L];

% シグマポイントに対応するyを計算
Y = mapcols(f, X')'; 

% yの期待値
ym = sum(W*Y)';

% 共分散行列
Yd  = bsxfun(@minus, Y, ym'); % 平均値の除去
Xd  = bsxfun(@minus, X, xm'); % 平均値の除去
Pyy = Yd'*W*Yd;
Pxy = Xd'*W*Yd;

end

