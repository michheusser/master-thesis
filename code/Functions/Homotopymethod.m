function ss = Homotopymethod (a,b,N,Ker,G) 
close all
clear all
clc


Ker = {'x^2 t*x', 't*x t^2'};
G = {'x-1.5*x^3', '1/6*x^2'};
a = 0;
b = 1;


syms x t p fg df ds1 dd
clc
tic
for i = 1:length(G)
    eval(['g' int2str(i), '= char(G(i))' ';'])
    d =eval (['g' int2str(i)])
    eval(['ff' int2str(i), '=d' ';'])
end
for j = 1 : N-1
    for i = 1 : length(G)
        ds1 = eval(['ff' int2str(i)]);
        df = ['f' int2str(i) int2str(j)]*p^j; eval (['ff' int2str(i), '=ds1+df' ';']);
    end
    for i = 1: length(G)
          
        kkerr = char(Ker(i));
        for ii = 1: length(G)
            ds2 = eval (['ff' int2str(ii)]);
            eval (['ff' int2str(ii) 't','=subs(ds2,x,t)' ';']);
            ds4 = eval (['ff' int2str(ii) 't']);
            ds3 = ['f' int2str(ii) 't'];
            kkerr = subs(kkerr,ds3,ds4);
        end
        ds1 = eval(['ff' int2str(i)]);
        gg = eval(['g' int2str(i)]);
        Homotopyequation = collect((ds1-gg)-p*int (kkerr,t,a,b),p);
        taylorhom1 = taylor(Homotopyequation,p,j+1)-taylor(Homotopyequation,p,j);
        H11 = subs(taylorhom1,{['f' int2str(i) int2str(j)],p},{0,1});
        eval(['f',int2str(i) int2str(j), '= -H11']);
        ds1 = eval(['ff' int2str(i)]);
        ds2 = ['f' int2str(i) int2str(j)];
        ds3 = eval (['f' int2str(i) int2str(j)]);
        eval (['ff' int2str(i) ' = collect(subs(ds1,ds2,ds3,x))' ';']);
    end
end
p=1;
for i=1:length(G)
    D = subs(eval(['ff' int2str(i)]),p,1);
    fprintf ('The approximate solution of %s is',['f' int2str(i) '(x)'])
    vpa(eval(D),10)
end
toc
end