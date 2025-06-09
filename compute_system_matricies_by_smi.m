%[�@�\] ������ԓ���@(Subspace Model Identification: SMI)�ɂ���ԕ������̍s��A�AB�AC�AD���v�Z����B
%[����]
% PUf, YfUfT, YfPiprepYfT�F�����������ōX�V�����s��
%�@�@�@�@�@�@�@�@�@�@�@ n�F�V�X�e������(���m�Ƃ���)
%�@�@�@�@�@�@�@�@�@�@�@ m�F���͐�
%�@�@�@�@�@�@�@�@�@�@�@ p�F�o�͐�
%�@�@�@�@�@�@�@�@�@�@�@ r�F�n���P���s��̍s�������߂邽�߂̃p�����[�^
%[�߂�l] ��ԕ������̍s��A�AB�AC�AD
function [A,B,C,D] = compute_system_matricies_by_smi(PUf, YfUfT, YfPiprepYfT, n, m, p, r)

    %���ْl����(Singular Value Decomposition: SVD)
    [EE,SS,~] = svd(YfPiprepYfT);
    
    En = EE(:,1:n);
    tempEn = En(1:p*(r-1),1:n);

    %�s��A�AC�̐���l
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
    
    %�s��B�CD�̐���l
    D = DB(1:p,:);
    B = DB(p+1:size(DB,1),:);
end