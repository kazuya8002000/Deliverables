%{
********************************************************************************************************************************************
�����@�F�K�������v�Z����֐�
���l�@�F�K�����̌v�Z���@�́u�V�X�e������̊�b�C�����d�@��w�o�ŋǁC2012�Cp�D28�v���Q�Ƃ����D
�쐬�ҁFR. Numata
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