function psi_n = normalize_psi(y_c,x_c,x_k,y_k)
psi_n=nan;
if (y_k-y_c)>0
    psi_n=atan((x_k-x_c)/(y_k-y_c));
end
if (x_k-x_c)>=0 & (y_k-y_c)<0
    psi_n=pi+atan(y/x);
end
if (x_k-x_c)<0 & (y_k-y_c)<0
    psi_n=-pi+atan(y/x);
end
if (x_k-x_c)>0 & (y_k-y_c)==0
    psi_n=pi/2;
end
if (x_k-x_c)<0 & (y_k-y_c)==0
    psi_n=-pi/2;
end
if psi_n<0
    psi_n=psi_n+2*pi;
end
end