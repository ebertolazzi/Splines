X = 0:0.1:1 ;
Y = 0:0.25:2 ;
[XX,YY] = ndgrid(X,Y) ;
ZZ = franke(XX,YY) ;

spline2d('pippo','bicubic',X,Y,ZZ) ;

surf(XX,YY,ZZ), view(145,-2), set(gca,'Fontsize',16) ;


X = 0:0.01:1 ;
Y = 0:0.025:2 ;
[XX,YY] = ndgrid(X,Y) ;

ZZ = zeros(size(XX)) ;
for i=1:size(XX,1)
	for j=1:size(XX,2)
		ZZ(i,j) = spline2d('pippo',[XX(i,j);YY(i,j)]) ;
	end
end

surf(XX,YY,ZZ), view(145,-2), set(gca,'Fontsize',16) ;
