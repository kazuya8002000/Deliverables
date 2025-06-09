%lq分解の関数
function [Ahat,Bhat,Chat,Dhat] = lqR4SID(Lf0,Lf2,Lf3,L22L22,n,m,p,r)

Imr       = eye(m*r,m*r);
Ip        = eye(p,p);
Zero_pr1p = zeros(p*(r-1),p);
maxindex = p * r - 1;
thresholdForSystemOrder  = 0.6;%0.5;%0.6, 0.2
%thresholdForSystemOrder2 = 0.6;%0.5;%0.6, 0.2
invLf0 = Lf0\Imr;
n_lq = n;

            %lqR4SIDAllTimesRecursions = toc(ctimeLqr4sid);
            
            [UU,SS,~] = svd(L22L22);

            %lqR4SIDAllTimesExObsMat = toc(ctimeLqr4sid);
            
            %システム次数を固有値から推定する．
%{          
            isCompletedLq = false;
            lqRatios = zeros(maxindex,1);
            n_lq = 0;
            for i = 1 : p*r
                index = i+1;

                if(index > p*r)
                    if(isCompletedLq == false)
                        n_lq = index;%i;
                    end
                    break;
                end

                ratio = SS(index,index)/SS(i,i);

                lqRatios(i,:) = ratio;
                
               if (ratio < thresholdForSystemOrder && isCompletedLq == false)
                   n_lq = index;%i;
                   isCompletedLq = true;
                  %break;
               end
            end
%}           
            %lqR4SIDEigenvalues(sidx : eidx,  k) = diag(SS);
            %lqR4SIDEigenRatios(sidx2: eidx2, k) = lqRatios;

            %lqR4SIDEstimatedSystemOrders(cnt, k) = n_lq;

            Ok = UU(:,1:n_lq);

            tempOk = Ok(1:p*(r-1),1:n_lq);
            Ahat   = pinv(tempOk)*Ok(p+1:p*r,1:n_lq);
            Chat   = Ok(1:p,1:n_lq);

            flowerL    = UU(:,n_lq+1:size(SS,1))';
            Z          = flowerL*Lf2*invLf0;
            [rowFL, ~] = size(flowerL);
            [rowZ , ~] = size(Z);

            LL = zeros(r * rowFL, r*p);
            XX = zeros(r *  rowZ, m  );
            for j=1:r
                tmp = j-1;
                LL(1 + rowFL* tmp: rowFL * j, :) = [flowerL(:, 1 + p * tmp: p * r) zeros(rowFL, p * tmp)];
                XX(1 + rowZ * tmp: rowZ  * j, :) =        Z(:, 1 + m * tmp: m * j);
            end

            DB = pinv(LL * [Ip       , zeros(p, n_lq); 
                            Zero_pr1p, tempOk ]) * XX;
            Dhat = DB(1:p,:);
            Bhat = DB(p+1:size(DB,1),:);

            Lp0 = Lf0;
            Lp2 = Lf2;
            Lp3 = Lf3;

            estimatedEigenvalue = eig(Ahat);

            realParts      = real(estimatedEigenvalue);
            imaginaryParts = imag(estimatedEigenvalue);

           % for j = 1:n_lq
           %    lqR4SIDAllRealParts(j,k)      = realParts(j,:);
           %    lqR4SIDAllImaginaryParts(j,k) = imaginaryParts(j,:);
           % end

            %estimatedSystem = ss(Ahat,Bhat,Chat,Dhat,[]);

            %estimatedOutput = lsim(estimatedSystem,testInput)';

            %trueOutput = lsim(trueSystem, testInput)';

            %averageOutput = mean(trueOutput, 2);

            %sumElement1 = zeros(p, initTimestep);
            %sumElement2 = zeros(p, initTimestep);
            %for l=1:initTimestep
            %    subtraction1 = trueOutput(:,l)-averageOutput;
            %    subtraction2 = estimatedOutput(:,l)-trueOutput(:,l);

            %   sumElement1(:,l) = subtraction1.^2;
            %    sumElement2(:,l) = subtraction2.^2;
            %end
end