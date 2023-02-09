function delta = ObjFunction(x)

%sample objective function


E =x(:,1);
P = 100;
L =1.5 ;
I =2e-7 ;

delta = (P.*L.^3)./(3*E.*I);
