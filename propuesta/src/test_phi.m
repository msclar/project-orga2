clear all;
[x,y]=meshgrid(0:0.1:4,0:0.1:4);

Phi=1500.*(x-2).*exp(-(x-2).^2.-(y-2).^2);

max(max(Phi))
min(min(Phi))

%surf(x,y,Phi);
pcolor(x,y,Phi);
shading interp;
