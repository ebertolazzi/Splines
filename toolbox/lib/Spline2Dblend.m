%>
%> MATLAB class wrapper for the underlying C++ class.
%>
%> The construction of the 2D spline is done as follows
%>
%> - instantiate the spline object
%>
%> ```{matlab}
%>   spl = Spline2Dblend( data );
%>   % or
%>   spl = Spline2Dblend();
%>   spl.build( data );
%> ```
%>
%> `data.type` is a string and can be any of
%>
%> |    kind      |   meaning                                                |
%> | :----------- | :------------------------------------------------------- |
%> | 'bilinear'   | piecewise linear in X e Y direction                      |
%> | 'cubic'      | piecewise cubic in X e Y direction                       |
%> | 'akima'      | piecewise cubic with Akima non oscillatory construction  |
%> | 'biquintic'  | piecewise quintic in X e Y direction                     |
%>
%> Evaluation is simple
%>
%> ```{matlab}
%>   Z = spl.eval( X, Y, s );
%> ```
%>
%> where `X` and `Y` are scalar or vector or matrix of the same size.
%> the result `Z` is scalar or vector or matrix of the same size of the inputs.
%>
%> **example of usage**
%>
%> ```{matlab}
%>   X = -2:0.01:2;
%>   Y = -2:0.01:2;
%>   [XX,YY] = ndgrid(X,Y);
%>   ZZ0  = peaks(XX,YY);
%>   ZZ1  = peaks(XX,YY);
%>   data.type = 'bicubic'
%>   data.surf0.xdata = XX
%>   data.surf0.ydata = YY
%>   data.surf0.zdata = ZZ0
%>   data.type = 'biquintic'
%>   data.surf1.xdata = XX
%>   data.surf1.ydata = YY
%>   data.surf1.zdata = ZZ0
%>   spl = Spline2Dblend(data);
%>
%>   surf(XX,YY,ZZ), view(145,-2), set(gca,'Fontsize',16);
%>
%>   X = -2:0.1:2;
%>   Y = -2:0.1:2;
%>   [XX,YY] = ndgrid(X,Y);
%>   s = 0.5;
%>
%>   ZZ = spl.eval(XX,YY,s);
%> ```
%>
%> @html_image{exampleSurf.png,width=80%}
%>
classdef Spline2Dblend < handle
  properties (SetAccess = private, Hidden = true)
    objectHandle; % Handle to the underlying C++ class instance
    call_delete;
  end

  methods
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Build a spline given a table of point and type:
    %>
    %>
    %> ```{matlab}
    %>   spl = Spline2Dblend( data ); % initialize
    %> ```
    %>
    function self = Spline2Dblend( data )
      self.objectHandle = Spline2DblendMexWrapper( 'new', data );
      self.call_delete  = true;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function delete(self)
      if self.call_delete
        Spline2DblendMexWrapper( 'delete', self.objectHandle );
      end
      self.call_delete  = false;
      self.objectHandle = 0;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Build a spline given a table of point and type:
    %>
    %>
    %> ```{matlab}
    %>   spl.build( data );
    %>   spl.build( 'file_name.json' );
    %> ```
    %>
    %> `X` and `Y` are vector `Z` is a matrix `nx x ny`
    %>
    %>
    function build( self, data )
      Spline2DblendMexWrapper( 'build', self.objectHandle, data );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Evaluate spline (and its derivatives) at (x,y)
    %>
    %>
    %> ```{matlab}
    %>   P = spl.eval( x, y, s );
    %>   [P,Px,Py] = spl.eval( x, y, s );
    %>   [P,Px,Py,Pxx,Pxy,Pyy] = spl.eval( x, y, s );
    %> ```
    %>
    %> `x` and `y` are salar or vector or matrix of the same size
    %>
    function varargout = eval( self, x, y, s )
      varargout{1} = Spline2DblendMexWrapper( 'eval', self.objectHandle, x, y, s );
      if nargout >= 3
        varargout{2} = Spline2DblendMexWrapper( 'eval_Dx', self.objectHandle, x, y, s );
        varargout{3} = Spline2DblendMexWrapper( 'eval_Dy', self.objectHandle, x, y, s );
      end
      if nargout >= 6
        varargout{4} = Spline2DblendMexWrapper( 'eval_Dxx', self.objectHandle, x, y, s );
        varargout{5} = Spline2DblendMexWrapper( 'eval_Dxy', self.objectHandle, x, y, s );
        varargout{6} = Spline2DblendMexWrapper( 'eval_Dyy', self.objectHandle, x, y, s );
      end
      if ~( nargout == 1 || nargout == 3 || nargout == 6 ) 
        error( 'Spline2Dblend.eval, nargout = %d must be 1, 3 or 6\n', nargout);
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Evaluate spline `x` derivative at (x,y)
    %>
    %>
    %> ```{matlab}
    %>   Dx = spl.eval_Dx( x, y, s );
    %> ```
    %>
    %> `x` and `y` are salar or vector or matrix of the same size
    %>
    function Dx = eval_Dx( self, x, y, s )
      Dx = Spline2DblendMexWrapper( 'eval_Dx', self.objectHandle, x, y, s );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Evaluate spline `y` derivative at (x,y)
    %>
    %>
    %> ```{matlab}
    %>   Dy = spl.eval_Dy( x, y, s );
    %> ```
    %>
    %> `x` and `y` are salar or vector or matrix of the same size
    %>
    function Dy = eval_Dy( self, x, y, s )
      Dy = Spline2DblendMexWrapper( 'eval_Dy', self.objectHandle, x, y, s );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Evaluate spline `x` second derivative at (x,y)
    %>
    %> ```{matlab}
    %>   Dxx = spl.eval_Dxx( x, y, s );
    %> ```
    %>
    %> `x` and `y` are salar or vector or matrix of the same size
    %>
    function Dxx = eval_Dxx( self, x, y, s )
      Dxx = Spline2DblendMexWrapper( 'eval_Dxx', self.objectHandle, x, y, s );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Evaluate spline `xy` second derivative at (x,y)
    %>
    %> ```{matlab}
    %>   Dxy = spl.eval_Dxy( x, y, s );
    %> ```
    %>
    %> `x` and `y` are salar or vector or matrix of the same size
    %>
    function Dxy = eval_Dxy( self, x, y, s )
      Dxy = Spline2DblendMexWrapper( 'eval_Dxy', self.objectHandle, x, y, s );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Evaluate spline `y` second derivative at (x,y)
    %>
    %> ```{matlab}
    %>   Dyy = spl.eval_Dyy( x, y, s );
    %> ```
    %>
    %> `x` and `y` are salar or vector or matrix of the same size
    %>
    function Dyy = eval_Dyy( self, x, y )
      Dyy = Spline2DblendMexWrapper( 'eval_Dyy', self.objectHandle, x, y, s );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Get definition domain of the 2D-spline
    %>
    %> ```{matlab}
    %>   res = spl.x_min();
    %> ```
    %>
    function res = x_min( self )
      res = Spline2DblendMexWrapper( 'x_min', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Get definition domain of the 2D-spline
    %>
    %> ```{matlab}
    %>   res = spl.x_max();
    %> ```
    %>
    function res = x_max( self )
      res = Spline2DblendMexWrapper( 'x_max', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Get definition domain of the 2D-spline
    %>
    %> ```{matlab}
    %>   res = spl.y_min();
    %> ```
    %>
    function res = y_min( self )
      res = Spline2DblendMexWrapper( 'y_min', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Get definition domain of the 2D-spline
    %>
    %> ```{matlab}
    %>   res = spl.y_max();
    %> ```
    %>
    function res = y_max( self )
      res = Spline2DblendMexWrapper( 'y_max', self.objectHandle );
    end
  end
end
