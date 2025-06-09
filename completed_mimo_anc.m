%�v���O�����̍ŏ��ł͊����̕ϐ��Ɨ�����K������������Ɨǂ��B
%�����̕�����ԓ���@
clc;   %�R�}���h�E�B���h�E�̗���������������B
clear; %�ϐ�������������B

%�V�~�����[�V��������
Fs                   = 1000;                         %�T���v�����O���g��
Ts                   = 1/Fs;                         %�T���v�����O����
max_time_s           = 60;                           %�V�~�����[�V��������(60�b = 1����)
block_length         = 1000;                         %�f�[�^�u���b�N�̒���
total_samples        = max_time_s*block_length;      %�f�[�^�̑���(60000)
time_axis            = (0:total_samples-1) / Fs;     %���Ԏ�
time_axis_J          = (0:max_time_s-1) ;            %�]���֐��̎��Ԏ�
num_blocks           = total_samples / block_length; %�f�[�^�u���b�N�̐�(max_time_s�Ɠ����͂�)
SNR_dB               = 36;                           %SN��
u_amp                = 100;                          %����u�̐U��
r                    = 5;%10;                           %������ԓ���@�̃n���P���s��̍s�������߂�p�����[�^
hankel_column_num    = block_length - r + 1;         %������ԓ���@�̃n���P���s��̗�
j_allowable_value_dB = -21;%-22;                     %�]���֐�J�̋��e�w�W
noise_variance       = 0.01;

%�v���C�}���p�X�̌W���s��̏����l
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

%�v���C�}���p�X�̃V�X�e������n�A���͐�m�A�o�͐�p
n = size(A_P0, 2);
m = size(B_P,  2);
p = size(C_P,  1);

%�ŗL�l�̏����l
[V0, D0] = eig(A_P0);

%�Z�J���_���p�X�̌W���s��
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

%�Z�J���_���p�X�̋t�V�X�e��
Ds_inv = Ds^-1;
As_inv = As - Bs*Ds_inv*Cs;
Bs_inv = Bs/Ds;
Cs_inv = -Ds_inv*Cs;

%�R���g���[���̌W���s��̏����l
A_C = zeros(5, 5);
B_C = zeros(5, 2);
C_C = zeros(2, 5);
D_C = zeros(2, 2);

%���肵���v���C�}���p�X�̌W���s��̏����l
est_A = zeros(3, 3);
est_B = zeros(3, 2);
est_C = zeros(2, 3);
est_D = zeros(2, 2);

%����u�Əo�͂̐^�ly(�G���������O�̏o��)
u =  u_amp * randn(m, total_samples);
%u = u_amp * idinput([total_samples, m], 'prbs')';
y = zeros(p, total_samples);
y_S = zeros(p, total_samples);
y_C = zeros(p, total_samples);
est_y_k = zeros(p, total_samples);
J = zeros(p, total_samples);
e_k_st = zeros(2, total_samples);
c_st = zeros(1, total_samples);

%�ϑ��G��
noise = sqrt(noise_variance).*randn(p, total_samples);

%��Ԃ̏����l
x_S_k = zeros(2,1);
x_P_k = zeros(n,1);
est_x_P_k = zeros(n,1);
inv_x_S_k = zeros(2,1);

%���肵���v���C�}���p�X��ۑ����邽�߂̃Z��
estimated_primary_paths = cell(1, num_blocks);

%�]���֐��̌v�Z�l��ۑ����邽�߂̔z��
j_values = zeros(1, max_time_s);

%�ȉ��̍s��Ƀv���C�}���p�X�̌ŗL�l��ۑ����āA�ǂ�ȓ���������V�X�e���Ȃ̂����Ă݂�B
eigs_primary = zeros(n, total_samples);

%1�b���Ƃɕ�����ԓ���@�����s����̂ŁA��O��for�����J�n1�A�I��60�Ƃ��Ă���B
for t_s = 1 : max_time_s
    
    e_k_st_J = zeros(2, block_length);
    g_J = zeros(p, block_length);
    
    %���U���Ԃ̃X�e�b�v��k���f�[�^�̃u���b�N����؂�ɂȂ�悤�ɁAfor���̍ŏ�(k_start)�ƍŌ�(k_end)���v�Z����B
    %(k_start�Ck_end) = (1�C1000)�A(1001�C2000)�A(2001�C3000)�A...�A(59001�C60000)�ƂȂ�悤�Ɍv�Z���Ă���B
    k_start = 1 + block_length * (t_s - 1);
    k_end   = block_length * t_s;

    %J�̃C���f�b�N�X��1000�ɌŒ肷�邽�߂̃C���f�b�N�X�̒�`
    index = 1;

    for k = k_start :  k_end
     
        %�v���C�}���p�X�̎��ςȌW���s��Ap(k)����邽�߂̑Ίp�s��
        tmp_coefficient = total_samples/2;
        Lpk = cos((pi/3) * (k - tmp_coefficient) / tmp_coefficient) * D0;
        c = cos((pi/3) * (k - tmp_coefficient) / tmp_coefficient);
        c_st(:, k) = c;
            
        %Ap(k)���v�Z����B�iApk = V0 * Lpk * inv(V0)�Ə����Ă��������ǁA���L�̏������̕����v�Z�������Đ��l���x�������B�j
        Apk = (V0 * Lpk) / V0;
            
        eigs_primary(:, k) = eig(Apk);
        
        %�v���C�}���p�X�̗��U���ԏ�ԕ�����(3)��
        x_P_k_1 = Apk * x_P_k + B_P * u(:, k);
        q_k     = C_P * x_P_k + D_P * u(:, k); %q(k)
        x_P_k   = x_P_k_1;
        y(:, k) = q_k + noise(:, k);
        
        %���肵���v���C�}���p�X�̗��U���ԏ�ԕ�����
        est_x_P_k_1 = est_A * est_x_P_k + est_B * u(:, k);
        est_y_P_k = est_C * est_x_P_k + est_D * u(:, k);
        
        est_y_k(:, k) = est_y_P_k;

        %�R���g���[���̗��U���ԏ�ԕ�����
        x_C = [est_x_P_k;
                  inv_x_S_k];%x_S_k];
        x_C_k_1 = A_C * x_C + B_C * u(:, k);
        y_C_k = C_C * x_C + D_C * u(:, k);
        y_C(:, k) = y_C_k;

        %���肵���v���C�}���p�X�̏�Ԃ̍X�V
        est_x_P_k = est_x_P_k_1;

        %�Z�J���_���p�X�̏�ԕ������̒�`
        x_S_k_1   = As * x_S_k + Bs * y_C_k;
        g_k       = Cs * x_S_k + Ds * y_C_k;%g(k)
        x_S_k     = x_S_k_1;
        y_S(:, k) = g_k;
        g_J(:, index) = g_k;
        
        %�Z�J���_���p�X�̋t�V�X�e���̏�ԕ�����
        inv_x_S_k_1 = As_inv * inv_x_S_k + Bs_inv * est_y_P_k;
        inv_y_S_k = Cs_inv * inv_x_S_k + Ds_inv * est_y_P_k;
        inv_x_S_k = inv_x_S_k_1;

        %�c������
        e_k = q_k - g_k;
        
        %�O���t�̃v���b�g�Ɏg���c�������̊i�[
        e_k_st(:, k) = e_k;

        %�]���p�����[�^J�Ɏg���c�������̊i�[
        e_k_st_J(:, index) = e_k;

        %�C���f�b�N�X�̍Ē�`
        index = index + 1;%�C���f�b�N�X��1001�ɂȂ�
         
    end

   J = 10*log10((sum(e_k_st_J(1, :).^2)+sum(e_k_st_J(2, :).^2))/(sum(g_J(1, :).^2)+sum(g_J(2, :).^2)));
   j_values(:, t_s) = J;

   if(J > j_allowable_value_dB)
%        for k = k_start :  k_end

        %�V�X�e������n�͊��m�Ƃ��āA������ԓ���@�����s���ăv���C�}���p�X�̏�ԕ������𐄒肷��B
       
        [est_A, est_B, est_C, est_D] = PerformSubspaceIdentificationByGivenSystemOrder(u(:, k_start:k_end+r-1), y(:, k_start:k_end+r-1), r, hankel_column_num, n);%hankel_column_num
        estimated_primary_paths{t_s} = ss(est_A, est_B, est_C, est_D, Ts);
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
        %�R���g���[���̌W���s����v�Z����B
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

%fit�̌v�Z
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

fit = ComputeFit(y(:, 1:60000), est_y_k);

disp(fit)

%�v���C�}���V�X�e���̓������O���t�Ƃ��ăv���b�g����B
fig_primary = figure('Color', 'White');
plot(time_axis, eigs_primary);
xlabel('Time (s)');
ylabel('Eigenvalues');
title('Eigenvalues of primary path');
%{
%���͂Əo�͂��O�̂��߃O���t������Ċm�F���Ă݂�B
fig_u = figure('Color', 'White');
plot(time_axis, u);
ylim([-1500, 1500]);
xlabel('Time (s)');
ylabel('Amplitude');
title('Input');

fig_y = figure('Color', 'White');
plot(time_axis, y);
ylim([-1500, 1500]);
xlabel('Time (s)');
ylabel('Amplitude');
title('Output');
%}
%�c���������O���t�Ƃ��ăv���b�g����B
fig_e_1 = figure('Color', 'White');
plot(time_axis,  y(1, :), "r");
ylim([-1500, 1500])
xlabel('Time (s)');
ylabel('Amplitude');
title('e_1');

hold on

plot(time_axis, e_k_st(1, :), "k");

hold off

fig_e_2 = figure('Color', 'White');
plot(time_axis,  y(2, :), "r");
ylim([-1500, 1500])
xlabel('Time (s)');
ylabel('Amplitude');
title('e_2');

hold on

plot(time_axis, e_k_st(2, :), "k");

hold off

%J���O���t�Ƃ��ăv���b�g����B
fig_J = figure('Color', 'White');
plot(time_axis_J, j_values);
ylim([-40, -10])
xlabel('Time (s)');
ylabel('J(dB)'); 
title('J');

%�Z�J���_���p�X�̏o�͂��O���t�Ƃ��ăv���b�g����B
fig_y_S = figure('Color', 'White');
plot(time_axis, y_S);
ylim([-1500, 1500])
xlabel('Time (s)');
ylabel('Eigenvalues');
title('y_S')
%{
fig_est_y_k = figure('Color', 'White');
plot(time_axis, est_y_k);
ylim([-1500, 1500])
xlabel('Time (s)');
ylabel('Eigenvalues');
title('est_y_k')

fig_y_C = figure('Color', 'White');
plot(time_axis, y_C);
ylim([-1500, 1500])
xlabel('Time (s)');
ylabel('Eigenvalues');
title('y_C')
%}
fig_cos = figure('Color', 'White');
plot(time_axis, c_st);
xlabel('Time (s)');
ylabel('cos');
xlim([0, 60])
title('cos');