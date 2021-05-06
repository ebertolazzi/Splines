%%-------------------------------------------------------------------------%%
%%                                                                         %%
%%  Copyright (C) 2020                                                     %%
%%                                                                         %%
%%         , __                 , __                                       %%
%%        /|/  \               /|/  \                                      %%
%%         | __/ _   ,_         | __/ _   ,_                               %%
%%         |   \|/  /  |  |   | |   \|/  /  |  |   |                       %%
%%         |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/                      %%
%%                           /|                   /|                       %%
%%                           \|                   \|                       %%
%%                                                                         %%
%%      Enrico Bertolazzi                                                  %%
%%      Dipartimento di Ingegneria Industriale                             %%
%%      Universita` degli Studi di Trento                                  %%
%%      email: enrico.bertolazzi@unitn.it                                  %%
%%                                                                         %%
%%-------------------------------------------------------------------------%%

close all;
clc;


x  = 1:1:100;
y  = x.^2;
yd = 2*x;
xx = 1:0.01:100;

fmt1 = {'LineWidth',3,'Color',[0.7,0.7,1]};
fmt2 = {'o','Color','red','MarkerSize',5,'MarkerFaceColor','yellow','MarkerEdgeColor','red'};
fmt3 = {'-','Color','black','Linewidth',2};
fmt4 = {'o','Color','blue','MarkerSize',3,'MarkerFaceColor','black','MarkerEdgeColor','black'};

figure();

plot( x, y, fmt2{:} ) ;
hold on;

S1 = Spline1D( 'cubic' );
S1.build( x, y );

yy = S1.eval( xx );
plot( xx, yy, fmt1{:} ) ;

S2 = Spline1D( 'hermite' );
S2.build( x, y, yd );
yy = S2.eval( xx );
plot( xx, yy, fmt3{:} ) ;
