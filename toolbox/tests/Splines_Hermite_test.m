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
%%      Università degli Studi di Trento                                   %%
%%      email: enrico.bertolazzi@unitn.it                                  %%
%%                                                                         %%
%%-------------------------------------------------------------------------%%

H = BaseHermite();

t = 0:1/1000:1;

[h0,h1,h2,h3] = H.base( t );
subplot(2,2,1);
plot( t, h0, t, h1, t, h2, t, h3 );

[h0,h1,h2,h3] = H.base_D( t );
subplot(2,2,2);
plot( t, h0, t, h1, t, h2, t, h3 );

[h0,h1,h2,h3] = H.base_DD( t );
subplot(2,2,3);
plot( t, h0, t, h1, t, h2, t, h3 );

[h0,h1,h2,h3] = H.base_DDD( t );
subplot(2,2,4);
plot( t, h0, t, h1, t, h2, t, h3 );
