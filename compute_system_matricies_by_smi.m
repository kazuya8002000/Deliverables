%[機能] 部分空間同定法(Subspace Model Identification: SMI)により状態方程式の行列A、B、C、Dを計算する。
%[引数]
% PUf, YfUfT, YfPiprepYfT：逐次方程式で更新した行列
%　　　　　　　　　　　 n：システム次数(既知とする)
%　　　　　　　　　　　 m：入力数
%　　　　　　　　　　　 p：出力数
%　　　　　　　　　　　 r：ハンケル行列の行数を決めるためのパラメータ
%[戻り値] 状態方程式の行列A、B、C、D
function [A,B,C,D] = compute_system_matricies_by_smi(PUf, YfUfT, YfPiprepYfT, n, m, p, r)

    %特異値分解(Singular Value Decomposition: SVD)
    [EE,SS,~] = svd(YfPiprepYfT);
    
    En = EE(:,1:n);
    tempEn = En(1:p*(r-1),1:n);

    %行列A、Cの推定値
    A = pinv(tempEn)*En(p+1:p*r,1:n);
    C = En(1:p,1:n);

    ALPHA = EE(:, n+1:size(SS,1))';
    BETA  = ALPHA*YfUfT*PUf;

    [rowAlpha, ~] = size(ALPHA);
    [rowBETA,  ~] = size(BETA);

    Tm_ALPHA = zeros(rowAlpha*r,p*r);
    BETA_T   = zeros(rowBETA*r,p);
    
    for l=1:r
        Tm_ALPHA(rowAlpha*(l-1)+1:rowAlpha*l,:) = [ALPHA(:,1+(l-1)*p:p*r) zeros(rowAlpha,p*(l-1))];
        BETA_T(rowBETA*(l-1)+1:rowBETA*l,:) = BETA(:,1+(l-1)*m:m*l);
    end

    B_En1 = [eye(p,p)        , zeros(p,n); 
             zeros(p*(r-1),p), tempEn   ];

    DB = pinv(Tm_ALPHA*B_En1)*BETA_T;
    
    %行列B，Dの推定値
    D = DB(1:p,:);
    B = DB(p+1:size(DB,1),:);
end