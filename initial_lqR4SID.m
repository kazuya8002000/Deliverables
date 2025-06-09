%初期値、大住
function [Lf0,Lf2,Lf3,L22L22] = initial_lqR4SID(uip,uif,yip,yif,Lp0,Lp2,Lp3,m,r)

Imr       = eye(m*r,m*r);

%逐次方程式
Lf0    = Lp0+uif*uif'-uip*uip';
invLf0 = Lf0\Imr;
Lf2    = Lp2+yif*uif'-yip*uip';
Lf3    = Lp3+yif*yif'-yip*yip';
L22L22 = Lf3-Lf2*invLf0*Lf2';
end