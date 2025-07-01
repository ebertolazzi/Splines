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

close all;

% Ottieni la dimensione dello schermo [left bottom width height]
scrn = get(0, 'ScreenSize');

% Calcola dimensioni all'80%
width  = round(scrn(3) * 0.8);
height = round(scrn(4) * 0.8);

% Centra la figura
left = round((scrn(3) - width) / 2);
bottom = round((scrn(4) - height) / 2);

% Crea la figura
figure('Position', [left bottom width height]);

V1 = zeros(10,11);
V1(3:6,4:7) = 3/7; % First step
V1(6:7,6:7) = [1/4 1/5; 1/5 1/4]; % Middle step
V1(8:9,8:9) = 1/2; % Last step
V1(1,1)     = -0.1;
V1(10,11)   = -0.1;
V1 = flip(V1,2);
X1 = 1:1:10;
Y1 = 1:1:11;


V2 = zeros(10,11);
V2(2:5,2:5) = 3/7; % First step
V2(6:7,6:7) = [1/4 1/5; 1/5 1/4]; % Middle step
V2(8:9,8:9) = 1/2; % Last step
V2(1,1)     = -0.1;
V2(10,11)   = -0.1;
V2 = flip(V2,2);
X2 = 1:1:10;
Y2 = 1:1:11;

type = { 'bilinear', 'bicubic', 'biquintic', 'akima' };
set(gca,'Fontsize',16);

data.surf0.spline_type     = 'bicubic';
data.surf0.xdata           = X1;
data.surf0.ydata           = Y1;
data.surf0.zdata           = V1;
data.surf0.fortran_storage = true;
data.surf0.transposed      = true;


data.surf1.spline_type     = 'biquintic';
data.surf1.xdata           = X2;
data.surf1.ydata           = Y2;
data.surf1.zdata           = V2;
data.surf1.fortran_storage = true;
data.surf1.transposed      = true;

S = Spline2Dblend(data);

x_min = S.x_min();
x_max = S.x_max();
y_min = S.y_min();
y_max = S.y_max();

X = x_min:0.01:x_max;
Y = y_min:0.01:y_max;
[XX,YY] = ndgrid(X,Y);


for k=1:4

  fprintf('Plot n.%d\n',k);

  ZZ = S.eval(XX,YY,(k-1)/3);

  if any(isnan(ZZ(:)))
    fprintf('Trovati NaN sulla superfice\n');
  end
  if any(isinf(ZZ(:)))
    fprintf('Trovati Inf sulla superfice\n');
  end

  subplot(2,2,k);

  surf(XX,YY,ZZ,'Linestyle',':');
  title(type{k});

  axis tight

  zlim([-0.5,0.6]);

  view(60,60);

end
