%初期値、奥
function [PUfEW,YfUfTEW,YfPiprepYfTEW] = initial_ewR4SID(PUpEW,YpUpTEW,YpPiperpYpTEW,uf,yf)

block_length         = 2000;%ここを変えた。長すぎるとよくない。
N0 = block_length;

%忘却係数
    forgettingFactor    = 1-1/N0;
    invForgettingFactor = 1/forgettingFactor;

%逐次方程式
    PUpuf         = PUpEW * uf;
    delta         = 1/(forgettingFactor + uf' * PUpuf);
    PUfEW         = invForgettingFactor * (PUpEW - (delta*PUpuf) * PUpuf');
    YfUfTEW       = forgettingFactor * YpUpTEW + yf * uf';
    q             = yf - YpUpTEW * PUpuf;
    YfPiprepYfTEW = forgettingFactor * (YpPiperpYpTEW + (delta * q) * q');
            