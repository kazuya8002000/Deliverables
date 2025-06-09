
clc;   %コマンドウィンドウの履歴を初期化する。
clear; %変数を初期化する。

%シミュレーション条件
Fs                   = 1000;                         %サンプリング周波数
Ts                   = 1/Fs;                         %サンプリング時間
max_time_s           = 60;                           %シミュレーション時間(60秒 = 1分間)
block_length         = 1000;                         %データブロックの長さ
total_samples        = max_time_s*block_length;      %データの総数(60000)
SNR_dB               = 36;                           %SN比
u_amp                = 100;                          %入力uの振幅
r                    = 5;                            %部分空間同定法のハンケル行列の行数を決めるパラメータ
hankel_column_num    = block_length - r + 1;         %部分空間同定法のハンケル行列の列数
initial_data_num     = block_length + r - 1;
noise_variance       = 0.01;
time_steps           = total_samples + initial_data_num;
time_axis            = (0:time_steps-1) / Fs;        %時間軸
time_axis_2          = (0:total_samples-1) / Fs;

%プライマリパスの係数行列の初期値
A_P0 = [ 0.5296, -0.4760,  0.1238;
        -0.4760, -0.0974,  0.1354;
         0.1238,  0.1354, -0.8233];

B_P = [-1.1470, -0.0376;
        1.1910,  0.3273;
             0,  0];

C_P = [-0.1867, -0.5883, -0.1364;
        0.7258,       0,  0.1139];

D_P = [1.0670, 0;
            0, 0.8];

system0 = ss(A_P0, B_P, C_P, D_P, []);

%プライマリパスのシステム次数n、入力数m、出力数p
n = size(A_P0, 2);
m = size(B_P,  2);
p = size(C_P,  1);

%入力uと出力の真値y(雑音が加わる前の出力)
u = u_amp * randn(m, time_steps);
%u = idinput([time_steps, m], 'prbs')';
y_P = zeros(p, time_steps);
y_P_plot = zeros(p, total_samples);
y_S = zeros(p, time_steps);
x_C_k = zeros(3, 1);
y_C = zeros(p, total_samples);
est_y_k = zeros(p, total_samples);
e_k_st = zeros(2, time_steps);
%Apk_st = cells();
c_st = zeros(1, time_steps);
svd_st = zeros(3, time_steps);
est_svd_st = zeros(3, time_steps);
norm_st = zeros(1, time_steps);
est_norm_st = zeros(1, time_steps);

%観測雑音
noise = sqrt(noise_variance).*randn(p, time_steps);

%プライマリパスの初期値の出力
u0 = u(:, 1:initial_data_num);
y0 = lsim(system0, u0)' + noise(1:initial_data_num);

%ハンケル型入出力データ行列
U0 = zeros(m*r,block_length);
Y0 = zeros(p*r,block_length);
for i = 1 : r
    U0(1+m*(i-1):m*i, :) = u(:,i:block_length+i-1);
    Y0(1+p*(i-1):p*i, :) = y0(:,i:block_length+i-1);
end

%逐次方程式の初期値を定義する。
PUp         = inv(U0*U0');
YpUpT       = Y0*U0';
YpPiperpYpT = Y0*(eye(block_length,block_length)-pinv(U0)*U0)*Y0';

%出力を格納するための配列
y = zeros(p, time_steps);
y(:, 1:initial_data_num) = y0;

%結果を格納するための配列
true_real_parts      = zeros(n, total_samples); %固有値の実部の真値
true_imag_parts      = zeros(n, total_samples); %固有値の虚部の真値
estimated_real_parts = zeros(n, total_samples); %固有値の実部の推定値
estimated_imag_parts = zeros(n, total_samples); %固有値の実部の推定値

%固有値の初期値
[V0, D0] = eig(A_P0);

%セカンダリパスの係数行列
a11 = 1;
a22 = 1;
a12 = 0.8;
a21 = 0.8;

As = zeros(2, 2);
Bs = eye(2, 2);
Cs = [  0, a12;
      a21,  0];
Ds = [a11, 0;
     0, a22];

%セカンダリパスの逆システム
Ds_inv = Ds^-1;
As_inv = As - Bs*Ds_inv*Cs;
Bs_inv = Bs/Ds;
Cs_inv = -Ds_inv*Cs;

%コントローラの係数行列の初期値
A_C = zeros(5, 5);
B_C = zeros(5, 2);
C_C = zeros(2, 5);
D_C = zeros(2, 2);

%状態の初期値
x_S_k = zeros(2,1);
x_P_k = zeros(n,1);
est_x_P_k = zeros(n,1);

%以下の行列にプライマリパスの固有値を保存して、どんな動きをするシステムなのか見てみる。
eigs_primary = zeros(n, total_samples);

%推定したプライマリパスの固有値を保存
est_eig = zeros(n, total_samples);

%プライマリパスの作成
for k = 1 : time_steps

    %現在の時刻を定義する。
    current_k = initial_data_num + k;

    %プライマリパスの時変な係数行列Ap(k)を作るための対角行列
    tmp_coefficient = time_steps/2;%120000/2
    Lpk = cos((pi/3) * (k - tmp_coefficient) / tmp_coefficient) * D0;
    c = cos((pi/3) * (k - tmp_coefficient) / tmp_coefficient);
    c_st(:, k) = c;

    %Ap(k)を計算する。（Apk = V0 * Lpk * inv(V0)と書いてもいいけど、下記の書き方の方が計算が早くて数値精度も高い。）
    Apk = (V0 * Lpk) / V0;

    eigs_primary(:, k) = eig(Apk);

    %プライマリパスの離散時間状態方程式(3)式
    x_P_k_1 = Apk * x_P_k + B_P * u(:, k);
    q_k     = C_P * x_P_k + D_P * u(:, k); %q(k)
    x_P_k   = x_P_k_1;
    y_P(:, k) = q_k + noise(:, k);
    %y_P_plot(:, k) = q_k;

    %推定値と比較するために、固有値の真値を配列へ保存しておく。
    true_eig = eig(Apk);
    true_real_parts(:, k) = real(true_eig);
    true_imag_parts(:, k) = imag(true_eig);
    
    [EE,SS,~] = svd(Apk);
    svd_st(:, k) = eig(SS);
    norm_st(:, k) = norm(Apk);

end

%入出力データベクトルを作成する。
uf = zeros(m*r, 1);
up = zeros(m*r, 1);
yf = zeros(p*r, 1);
yp = zeros(p*r, 1);

for k = 1 : total_samples
 
    %現在の時刻を定義する。
    current_k = initial_data_num + k;

    y(:, current_k) = y_P(:, current_k);

    for i=1:r
        ufTemp            = u(:, block_length+k:current_k);
        upTemp            = u(:, k:r+(k-1));
        uf(m*i-1: m*i, :) = ufTemp(:, i);
        up(m*i-1: m*i, :) = upTemp(:, i);

        yfTemp            = y(:, block_length+k:current_k);
        ypTemp            = y(:, k:r+(k-1));
        yf(p*i-1: p*i, :) = yfTemp(:, i);
        yp(p*i-1: p*i, :) = ypTemp(:, i);
    end

    %逐次部分空間同定法を実行する。
    [PUfEW, YfUfTEW, YfPiprepYfTEW] = initial_ewR4SID(PUp, YpUpT, YpPiperpYpT, uf, yf);
    [est_A,est_B,est_C,est_D] = ewR4SID(PUfEW, YfUfTEW, YfPiprepYfTEW, n, m, p, r);

    %更新前の行列を再定義する。
    PUp         = PUfEW;
    YpPiperpYpT = YfPiprepYfTEW;
    YpUpT       = YfUfTEW;

    %推定値を保存する。
    est_eig(:, k) = eig(est_A);
    est_eig_k = eig(est_A);
    estimated_real_parts(:, k) = real(est_eig_k);
    estimated_imag_parts(:, k) = imag(est_eig_k);
    [EE,SS,~] = svd(est_A);
    est_svd_st(:, k) = eig(SS);
    est_norm_st(:, k) = norm(est_A);

    %推定したプライマリパスの離散時間状態方程式
    est_x_P_k_1 = est_A * est_x_P_k + est_B * u(:, k);
    est_y_P_k = est_C * est_x_P_k + est_D * u(:, k);
    
    est_y_k(:, k) = est_y_P_k;

    %コントローラの係数行列を計算する。
    element_A_C = Bs_inv * est_C;
    A_C = [est_A, zeros(3, 2);
         element_A_C, As_inv];
    B_C = [est_B;
        Bs_inv * est_D];
    C_C = [Ds_inv * est_C, Cs_inv];
    D_C = Ds_inv * est_D;

    %コントローラの離散時間状態方程式
    x_C = [est_x_P_k;
              x_S_k];%x_inv_S_k];
    x_C_k_1 = A_C * x_C + B_C * u(:, k);
    y_C_k = C_C * x_C + D_C * u(:, k);
    y_C(:, k) = y_C_k;

    %推定したプライマリパスの状態の更新
    est_x_P_k = est_x_P_k_1;
    
    %セカンダリパスの状態方程式の定義
    x_S_k_1   = As * x_S_k + Bs * y_C_k;
    g_k       = Cs * x_S_k + Ds * y_C_k;%g(k)
    x_S_k     = x_S_k_1;
    y_S(:, k) = g_k;
    
    %残留騒音
    e_k = y(:, k) - g_k;
    
    %グラフのプロットに使う残留騒音の格納
    e_k_st(:, k) = e_k;

end

%fitの表示

fit = ComputeFit(y(:, 1:60000), est_y_k);

disp(fit)

%残留騒音をグラフとしてプロットする。
fig_e_1 = figure('Color', 'White');
plot(time_axis,  y_P(1, :), "r");
ylim([-1500, 1500])
xlim([0, 60])
xlabel('Time (s)');
ylabel('Amplitude');
title('e_1');

hold on

plot(time_axis, e_k_st(1, :), "k");

hold off

fig_e_2 = figure('Color', 'White');
plot(time_axis,  y_P(2, :), "r");
ylim([-1500, 1500])
xlim([0, 60])
xlabel('Time (s)');
ylabel('Amplitude');
title('e_2');

hold on

plot(time_axis, e_k_st(2, :), "k");

hold off

%プライマリシステムの動きをグラフとしてプロットする。
fig_primary = figure('Color', 'White');
plot(time_axis, eigs_primary);
xlabel('Time (s)');
ylabel('Eigenvalues');
xlim([0, 60])
title('Eigenvalues of primary path');

%推定したプライマリシステムの動きをグラフとしてプロットする。
fig_estprimary1 = figure('Color', 'White');
plot(time_axis_2, est_eig);
xlabel('Time (s)');
ylabel('Eigenvalues');
xlim([0, 60])
title('Eigenvalues of estimated primary path');

%推定したプライマリシステムの動きをグラフとしてプロットする。
fig_estprimary2 = figure('Color', 'White');
plot(time_axis_2, estimated_real_parts);
xlabel('Time (s)');
ylabel('Eigenvalues');
ylim([-1, 0.8])
title('Eigenvalues of estimated primary path real');

fig_estprimary3 = figure('Color', 'White');
plot(time_axis_2, estimated_imag_parts);
xlabel('Time (s)');
ylabel('Eigenvalues');
ylim([-1, 0.8])
title('Eigenvalues of estimated primary path imag');
%{
g = get(gca);
g.Children

g = [(e_k_st(2, :));
      (y_P(2, :))];

set(gca, 'Children', [g.Children(2) g.Children(1)])
%}
fig_cos = figure('Color', 'White');
plot(time_axis, c_st);
xlabel('Time (s)');
ylabel('cos');
xlim([0, 60])
title('cos');

fig_svd = figure('Color', 'White');
plot(time_axis, svd_st);
xlabel('Time (s)');
ylabel('svd');
xlim([0, 60])
title('singular value decomposition of primarypath');

fig_est_svd = figure('Color', 'White');
plot(time_axis, est_svd_st);
xlabel('Time (s)');
ylabel('est_svd');
xlim([0, 60])
ylim([0, 1.25])
title('estimated singular value decomposition of secondarypath');

fig_norm = figure('Color', 'White');
plot(time_axis, norm_st);
xlabel('Time (s)');
ylabel('norm');
xlim([0, 60])
ylim([0, 1.3])
title('norm of primarypath');

fig_est_norm = figure('Color', 'White');
plot(time_axis, est_norm_st);
xlabel('Time (s)');
ylabel('estimated norm');
xlim([0, 60])
ylim([0, 1.3])
title('estimated norm of secondarypath');