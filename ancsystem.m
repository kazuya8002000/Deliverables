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
j_allowable_value_dB = -21;%-22;                     %評価関数Jの許容指標
noise_variance       = 0.01;
time_steps           = total_samples + initial_data_num;
time_axis            = (0:time_steps-1) / Fs;        %時間軸
time_axis_2          = (0:total_samples-1) / Fs;
time_axis_J          = (0:max_time_s) ;            %評価関数の時間軸(ここが元とちょっと違う)

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
inv_x_S_k = zeros(2,1);

%評価関数の計算値を保存するための配列
j_values = zeros(1, max_time_s);

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

%ハンケル型入出力データ行列のLQ分解
[~,R0] = qr([U0;Y0]',0);
L0     = R0';
L11_0  = L0(1:m*r,1:m*r);
L21_0  = L0(m*r+1:m*r+p*r,1:m*r);
L22_0  = L0(m*r+1:m*r+p*r,m*r+1:m*r+p*r);

%逐次方程式の初期値
Lp0 = L11_0*L11_0';
Lp2 = L21_0*L11_0';
Lp3 = L21_0*L21_0'+L22_0*L22_0';

%逐次方程式の初期値を定義する。
PUp         = inv(U0*U0');
YpUpT       = Y0*U0';
YpPiperpYpT = Y0*(eye(block_length,block_length)-pinv(U0)*U0)*Y0';

PUpEW         = inv(U0*U0');
YpUpTEW       = Y0*U0';
YpPiperpYpTEW = Y0*(eye(block_length,block_length)-pinv(U0)*U0)*Y0';

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

%推定したプライマリパスの係数行列の初期値
est_A = zeros(3, 3);
est_B = zeros(3, 2);
est_C = zeros(2, 3);
est_D = zeros(2, 2);

%状態の初期値
x_S_k = zeros(2,1);
x_P_k = zeros(n,1);
est_x_P_k = zeros(n,1);

%以下の行列にプライマリパスの固有値を保存して、どんな動きをするシステムなのか見てみる。
eigs_primary = zeros(n, total_samples);

%推定したプライマリパスの固有値を保存
est_eig = zeros(n, total_samples);

%1秒ごとに部分空間同定法を実行するので、大外のfor文を開始1、終了60としている。
for t_s = 1 : max_time_s + 1
    
    e_k_st_J = zeros(2, block_length);
    g_J = zeros(p, block_length);
    
    %離散時間のステップ数kがデータのブロック長区切りになるように、for分の最初(k_start)と最後(k_end)を計算する。
    %(k_start，k_end) = (1，1000)、(1001，2000)、(2001，3000)、...、(59001，60000)となるように計算している。
    k_start = 1 + block_length * (t_s - 1);
    k_end   = block_length * t_s;

    %Jのインデックスを1000個に固定するためのインデックスの定義
    index = 1;

    for k = k_start :  k_end
     
        %プライマリパスの時変な係数行列Ap(k)を作るための対角行列
        tmp_coefficient = total_samples/2;
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
        
        %推定したプライマリパスの離散時間状態方程式
        est_x_P_k_1 = est_A * est_x_P_k + est_B * u(:, k);
        est_y_P_k = est_C * est_x_P_k + est_D * u(:, k);
        
        est_y_k(:, k) = est_y_P_k;

        %コントローラの離散時間状態方程式
        x_C = [est_x_P_k;
                  inv_x_S_k];%x_S_k];
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
        g_J(:, index) = g_k;
        
        %セカンダリパスの逆システムの状態方程式
        inv_x_S_k_1 = As_inv * inv_x_S_k + Bs_inv * est_y_P_k;
        inv_y_S_k = Cs_inv * inv_x_S_k + Ds_inv * est_y_P_k;
        inv_x_S_k = inv_x_S_k_1;

        %残留騒音
        e_k = q_k - g_k;
        
        %グラフのプロットに使う残留騒音の格納
        e_k_st(:, k) = e_k;

        %評価パラメータJに使う残留騒音の格納
        e_k_st_J(:, index) = e_k;

        %インデックスの再定義
        index = index + 1;%インデックスが1001になる
         
    end

   J = 10*log10((sum(e_k_st_J(1, :).^2)+sum(e_k_st_J(2, :).^2))/(sum(g_J(1, :).^2)+sum(g_J(2, :).^2)));
   j_values(:, t_s) = J;

   if(J > j_allowable_value_dB)
%        for k = k_start :  k_end

        %システム次数nは既知として、部分空間同定法を実行してプライマリパスの状態方程式を推定する。
       
        [est_A, est_B, est_C, est_D] = PerformSubspaceIdentificationByGivenSystemOrder(u(:, k_start:k_end), y_P(:, k_start:k_end), r, hankel_column_num, n);
        %estimated_primary_paths{t_s} = ss(est_A, est_B, est_C, est_D, Ts);
%{ 
        sys = n4sid(u(:, k_start:k_end), y(:, k_start:k_end));
        est_A = sys.A;
        est_B = sys.B;
        est_C = sys.C;
        est_D = sys.D;

        [Result] = DSSI(y(:, k_start:k_end), u(:, k_start:k_end), Ts, 500);
        Result.Matrices.A = est_A;
        Result.Matrices.B = est_B;
        Result.Matrices.C = est_C;
        Result.Matrices.D = est_D;

        cut = 4;
        Dt=1/Ts;
        DAT = iddata(y(:, k)', u(:, k)', Dt);     % Input Output data
        TH1 = n4sid(DAT,cut);%DSSI
        est_A = TH1.A;                           % A matrix
        est_B = TH1.B;                           % B matrix
        est_C = TH1.C;                           % C matrix
        est_D = TH1.D;                           % D matrix
%}
        %コントローラの係数行列を計算する。
        element_A_C = Bs_inv * est_C;
        A_C = [est_A, zeros(3, 2);
             element_A_C, As_inv];
        B_C = [est_B;
            Bs_inv * est_D];
        C_C = [Ds_inv * est_C, Cs_inv];
        D_C = Ds_inv * est_D;

 %       end
    end
end

%fitの計算
%{
for k = 1 : total_samples

    ave_y_P_1(:, k) = sum(y(1, :))/total_samples;
    ave_y_P_2(:, k) = sum(y(2, :))/total_samples;
    %ave_y_S = sum(y_S)/total_samples;
    fit_y_1(:, k) = ((ave_y_P_1(:, k) - y(1, k))^2);
    fit_y_2(:, k) = ((ave_y_P_2(:, k) - y(2, k))^2);
    fit_y_C_1(:, k) = ((ave_y_P_1(:, k) - y_S(1, k))^2);
    fit_y_C_2(:, k) = ((ave_y_P_2(:, k) - y_S(2, k))^2);
end

fix_function_1 = (1-sqrt(sum(fit_y_C_1)/sum(fit_y_1))) * 100;
fix_function_2 = (1-sqrt(sum(fit_y_C_2)/sum(fit_y_2))) * 100;
%}

fit = ComputeFit(y_P(:, 1:60000), est_y_k);

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

%Jをグラフとしてプロットする。
fig_J = figure('Color', 'White');
plot(time_axis_J, j_values);
ylim([-40, -10])
xlabel('Time (s)');
ylabel('J(dB)'); 
title('J');

%初期値
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
inv_x_S_k = zeros(2,1);
g_k = zeros(2,1);

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
est_x_P_k_1 = zeros(3,1);
est_y_P_k = zeros(2,1);
x_C_k_1 = zeros(5,1);
y_C_k = zeros(2,1);
x_S_k_1 = zeros(2,1);
inv_x_S_k_1 = zeros(2,1);

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
    [PUfEW, YfUfTEW, YfPiprepYfTEW] = initial_ewR4SID(PUpEW, YpUpTEW, YpPiperpYpTEW, uf, yf);
    [est_A,est_B,est_C,est_D] = ewR4SID(PUfEW, YfUfTEW, YfPiprepYfTEW, n, m, p, r);

    %更新前の行列を再定義する。
    PUpEW         = PUfEW;
    YpPiperpYpTEW = YfPiprepYfTEW;
    YpUpTEW       = YfUfTEW;

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

%初期値
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
inv_x_S_k = zeros(2,1);
g_k = zeros(2,1);

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
est_x_P_k_1 = zeros(3,1);
est_y_P_k = zeros(2,1);
x_C_k_1 = zeros(5,1);
y_C_k = zeros(2,1);
x_S_k_1 = zeros(2,1);
inv_x_S_k_1 = zeros(2,1);

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
    %[PUf, YfUfT, YfPiprepYfT] = single_mil_based_r4sid(PUp, YpUpT, YpPiperpYpT, up, yp, uf, yf);
    %[est_A,est_B,est_C,est_D] = compute_system_matricies_by_smi(PUf, YfUfT, YfPiprepYfT, n, m, p, r);
    [Lf0,Lf2,Lf3,L22L22] = initial_lqR4SID(up,uf,yp,yf,Lp0,Lp2,Lp3,m,r);
    [est_A,est_B,est_C,est_D] = lqR4SID(Lf0,Lf2,Lf3,L22L22,n,m,p,r);

    %更新前の行列を再定義する。
    Lp0 = Lf0;
    Lp2 = Lf2;
    Lp3 = Lf3;

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

%初期値
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
inv_x_S_k = zeros(2,1);
g_k = zeros(2,1);

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
est_x_P_k_1 = zeros(3,1);
est_y_P_k = zeros(2,1);
x_C_k_1 = zeros(5,1);
y_C_k = zeros(2,1);
x_S_k_1 = zeros(2,1);
inv_x_S_k_1 = zeros(2,1);

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

    %逐次部分空間同定法を実行する。沼田先輩の！！single...で逐次方程式で使うやつを出して、compute,,,でABCDを出す。
    [PUf, YfUfT, YfPiprepYfT] = single_mil_based_r4sid(PUp, YpUpT, YpPiperpYpT, up, yp, uf, yf);
    [est_A,est_B,est_C,est_D] = compute_system_matricies_by_smi(PUf, YfUfT, YfPiprepYfT, n, m, p, r);
    
    %更新前の行列を再定義する。
    PUp         = PUf;
    YpUpT       = YfUfT;
    YpPiperpYpT = YfPiprepYfT;
    
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

%初期値の定義
   
    P_p     = eye(m*(n+p),m*(n+p));
    Ptp     = eye(n,n);

[theta_p,varphi2_p,Wtp,Lp,rotationCount,est_A,est_B,est_C,est_D] = InitialValue_Recursive_subspace_identification(U0,Y0,u,n,m,p,r);

%初期値
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
inv_x_S_k = zeros(2,1);
g_k = zeros(2,1);

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
est_x_P_k_1 = zeros(3,1);
est_y_P_k = zeros(2,1);
x_C_k_1 = zeros(5,1);
y_C_k = zeros(2,1);
x_S_k_1 = zeros(2,1);
inv_x_S_k_1 = zeros(2,1);

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

    %varphi2 = varphi2_p*kron(eye(m,m),A_past)+kron(u(:,current_k-1)',C_past);
    %varphi  = [kron(u(:,current_k)',eye(p,p)) varphi2];

    %逐次部分空間同定法を実行する。
    
    %system0 = ss(A_past0,B_past0,C_past0,D_past0,[]);
    %est_y0 = lsim(system0,u)';
    
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

    u_tver_1 = u(:,current_k-1);
    u_tver = u(:,current_k);
    ynew = y_P(:,current_k);%ちがうかも
    
    [Wtf,H,Ptf] = Compute_past(uf,yf,Wtp,Lp,Ptp,n,m,p,r);
    [est_A,est_B,est_C,est_D,P_f,theta,varphi2] = Recursive_subspace_identification(Wtf,theta_p,varphi2_p,P_p,u_tver,u_tver_1,ynew,n,m,p,r);
    
    theta_p   = theta;
    varphi2_p = varphi2;
    P_p       = P_f;
    Ptp       = Ptf;
    Wtp       = Wtf;
    Lp        = H(:,1:m*r+p*r);
    
    %推定値を保存する。
    %est_eig(:, k) = eig(est_A);
    %est_eig_k = eig(est_A);
    %estimated_real_parts(:, k) = real(est_eig_k);
    %estimated_imag_parts(:, k) = imag(est_eig_k);
    %[EE,SS,~] = svd(est_A);
    %est_svd_st(:, k) = eig(SS);
    %est_norm_st(:, k) = norm(est_A);

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