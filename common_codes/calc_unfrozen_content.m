%mv: total soil water content in m3/m3
%T: soil temperature in degC

function mvu=calc_unfrozen_content(mv,T)


C=1;
A_C=0.00306;
A_mv=0.394;
B_C=0.00582;
B_mv=-1.073;


A = A_C*C + A_mv*mv;
B = B_C*C + B_mv*mv;
mvu=A.*abs(T).^B;


end