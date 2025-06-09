%{
********************************************************************************************************************************************
説明　：適合率を計算する関数
備考　：適合率の計算方法は「システム同定の基礎，東京電機大学出版局，2012，p．28」を参照した．
作成者：R. Numata
********************************************************************************************************************************************
%}

function fit = ComputeFit(trueOutput, estimatedOutput)
    
    averageOutput = mean(trueOutput, 2);
    [outputCount, dataCount] = size(trueOutput);
    
    denoms = zeros(outputCount, dataCount);
    numes  = zeros(outputCount, dataCount);
    for columnNumber = 1 : dataCount
        denom = trueOutput(:, columnNumber) - averageOutput;
        nume  = estimatedOutput(:,columnNumber) - trueOutput(:,columnNumber);
        denoms(:, columnNumber) = denom.^2;
        numes(:, columnNumber)  = nume.^2;
    end
    denominator=sqrt(sum(denoms, 2));
    numerator=sqrt(sum(numes, 2));
    
    fit = (ones(outputCount, 1) - numerator ./ denominator) * 100;
end