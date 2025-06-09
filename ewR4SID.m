%奥先生
function [Ahat,Bhat,Chat,Dhat] = ewR4SID(PUfEW, YfUfTEW, YfPiprepYfTEW, n, m, p, r)

%Imr       = eye(m*r,m*r);
Ip        = eye(p,p);
Zero_pr1p = zeros(p*(r-1),p);
%maxindex = p * r - 1;
%thresholdForSystemOrder  = 0.6;%0.5;%0.6, 0.2
%thresholdForSystemOrder2 = 0.6;%0.5;%0.6, 0.2
%invLf0 = Lf0\Imr;
n_ew = n;

%ewR4SIDAllTimesRecursions = toc(ctimeEwr4sid);

[EE,SS,~] = svd(YfPiprepYfTEW);

%ewR4SIDAllTimesExObsMat = toc(ctimeEwr4sid);

%システム次数を固有値から推定する．
%isCompletedEw = false;
%ewRatios = zeros(maxindex,1);
%{            
n_ew = 0;
for i = 1 : p*r
    index = i+1;

    if(index > p*r)
        if(isCompletedEw == false)
            n_ew = index;%i;
        end
        break;
    end

    ratio = SS(index,index)/SS(i,i);

    %ewRatios(i,:) = ratio;
    
   if (ratio < thresholdForSystemOrder2 && isCompletedEw == false)
       n_ew = index;%i;
       isCompletedEw = true;
      %break;
   end
end
%}
%ewR4SIDEigenvalues(sidx  : eidx , k) = diag(SS);
%ewR4SIDEigenRatios(sidx2 : eidx2, k) = ewRatios;

%ewR4SIDEstimatedSystemOrders(cnt,k) = n_ew;

En = EE(:,1:n_ew);

tempEn = En(1:p*(r-1),1:n_ew);

Ahat = pinv(tempEn)*En(p+1:p*r,1:n_ew);
Chat = En(1:p,1:n_ew);

ALPHA = EE(:, n_ew+1:size(SS,1))';
BETA  = ALPHA*YfUfTEW*PUfEW;

[rowAlpha, columnAlpha] = size(ALPHA);
[rowBETA,  columnBETA]  = size(BETA);

Tm_ALPHA = zeros(rowAlpha*r,p*r);
BETA_T   = zeros(rowBETA*r,p);
for l=1:r
    Tm_ALPHA(rowAlpha*(l-1)+1:rowAlpha*l,:) = [ALPHA(:,1+(l-1)*p:p*r) zeros(rowAlpha,p*(l-1))];
    BETA_T(rowBETA*(l-1)+1:rowBETA*l,:) = BETA(:,1+(l-1)*m:m*l);
end

B_En1 = [       Ip, zeros(p,n_ew); 
         Zero_pr1p, tempEn        ];

DB   = pinv(Tm_ALPHA*B_En1)*BETA_T;
Dhat = DB(1:p,:);
Bhat = DB(p+1:size(DB,1),:);

%PUpEW         = PUfEW;
%YpPiperpYpTEW = YfPiprepYfTEW;
%YpUpTEW       = YfUfTEW;
