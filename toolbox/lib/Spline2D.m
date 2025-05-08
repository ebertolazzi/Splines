%>
%> MATLAB class wrapper for the underlying C++ class.
%>
%> The construction of the 2D spline is done as follows
%>
%> - instantiate the spline object
%>
%> ```{matlab}
%>   spl = Spline2D( kind, X, Y, ZZ );
%>   % or
%>   spl = Spline2D( kind );
%>   spl.build( X, Y, ZZ );
%> ```
%>
%> `kind` is a string and can be any of 
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
%>   Z = spl.eval( X, Y );
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
%>   ZZ  = peaks(XX,YY);
%>   spl = Spline2D('bicubic',X,Y,ZZ);
%>
%>   surf(XX,YY,ZZ), view(145,-2), set(gca,'Fontsize',16);
%>
%>   X = -2:0.1:2;
%>   Y = -2:0.1:2;
%>   [XX,YY] = ndgrid(X,Y);
%>
%>   ZZ = spl.eval(XX,YY);
%> ```
%>
%> @html_image{exampleSurf.png,width=80%}
%>
classdef Spline2D < handle
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
    %>   spl = Spline2D( kind );          % initialize
    %>   spl = Spline2D( kind, X, Y, Z ); % initialize and build
    %>   spl = Spline2D( file );          % reading data from file (yaml,json or toml)
    %> ```
    %>
    function self = Spline2D( name, varargin )
      self.objectHandle = Spline2DMexWrapper( 'new', name );
      self.call_delete  = true;
      if nargin > 1
        Spline2DMexWrapper( 'build', self.objectHandle, varargin{:} );
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function delete(self)
      if self.call_delete
        Spline2DMexWrapper( 'delete', self.objectHandle );
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
    %>   spl.build( X, Y, Z );
    %>   spl.build( 'file_name.json' );
    %> ```
    %>
    %> `X` and `Y` are vector `Z` is a matrix `nx x ny`
    %>
    %>
    function build( self, varargin )
      Spline2DMexWrapper( 'build', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Evaluate spline (and its derivatives) at (x,y)
    %>
    %>
    %> ```{matlab}
    %>   P = spl.eval( x, y );
    %>   [P,Px,Py] = spl.eval( x, y );
    %>   [P,Px,Py,Pxx,Pxy,Pyy] = spl.eval( x, y );
    %> ```
    %>
    %> `x` and `y` are salar or vector or matrix of the same size
    %>
    function varargout = eval( self, x, y )
      varargout{1} = Spline2DMexWrapper( 'eval', self.objectHandle, x, y );
      if nargout >= 3
        varargout{2} = Spline2DMexWrapper( 'eval_Dx', self.objectHandle, x, y );
        varargout{3} = Spline2DMexWrapper( 'eval_Dy', self.objectHandle, x, y );
      end
      if nargout >= 6
        varargout{4} = Spline2DMexWrapper( 'eval_Dxx', self.objectHandle, x, y );
        varargout{5} = Spline2DMexWrapper( 'eval_Dxy', self.objectHandle, x, y );
        varargout{6} = Spline2DMexWrapper( 'eval_Dyy', self.objectHandle, x, y );
      end
      if ~( nargout == 1 || nargout == 3 || nargout == 6 ) 
        error( 'Spline2D.eval, nargout = %d must be 1, 3 or 6\n', nargout);
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Evaluate spline `x` derivative at (x,y)
    %>
    %>
    %> ```{matlab}
    %>   Dx = spl.eval_Dx( x, y );
    %> ```
    %>
    %> `x` and `y` are salar or vector or matrix of the same size
    %>
    function Dx = eval_Dx( self, x, y )
      Dx = Spline2DMexWrapper( 'eval_Dx', self.objectHandle, x, y );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Evaluate spline `y` derivative at (x,y)
    %>
    %>
    %> ```{matlab}
    %>   Dy = spl.eval_Dy( x, y );
    %> ```
    %>
    %> `x` and `y` are salar or vector or matrix of the same size
    %>
    function Dy = eval_Dy( self, x, y )
      Dy = Spline2DMexWrapper( 'eval_Dy', self.objectHandle, x, y );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Evaluate spline `x` second derivative at (x,y)
    %>
    %> ```{matlab}
    %>   Dxx = spl.eval_Dxx( x, y );
    %> ```
    %>
    %> `x` and `y` are salar or vector or matrix of the same size
    %>
    function Dxx = eval_Dxx( self, x, y )
      Dxx = Spline2DMexWrapper( 'eval_Dxx', self.objectHandle, x, y );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Evaluate spline `xy` second derivative at (x,y)
    %>
    %> ```{matlab}
    %>   Dxy = spl.eval_Dxy( x, y );
    %> ```
    %>
    %> `x` and `y` are salar or vector or matrix of the same size
    %>
    function Dxy = eval_Dxy( self, x, y )
      Dxy = Spline2DMexWrapper( 'eval_Dxy', self.objectHandle, x, y );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Evaluate spline `y` second derivative at (x,y)
    %>
    %> ```{matlab}
    %>   Dyy = spl.eval_Dyy( x, y );
    %> ```
    %>
    %> `x` and `y` are salar or vector or matrix of the same size
    %>
    function Dyy = eval_Dyy( self, x, y )
      Dyy = Spline2DMexWrapper( 'eval_Dyy', self.objectHandle, x, y );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Set spline surface as closed in `x direction
    %>
    %> ```{matlab}
    %>   spl.make_x_closed();
    %> ```
    %>
    function make_x_closed( self )
      Spline2DMexWrapper( 'make_x_closed', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Set spline surface as opened in `x direction
    %>
    %> ```{matlab}
    %>   spl.make_x_opened();
    %> ```
    %>
    function make_x_opened( self )
      Spline2DMexWrapper( 'make_x_opened', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Return true if spline surface is of closed type in `x`-direction
    %>
    %> ```{matlab}
    %>   ok = spl.is_x_closed();
    %> ```
    %>
    function ok = is_x_closed( self )
      ok = Spline2DMexWrapper( 'is_x_closed', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Make spline surface computable only in the `x`-range where is defined.
    %> If `x` is outside range an error is produced.
    %>
    %> ```{matlab}
    %>   spl.make_x_bounded();
    %> ```
    %>
    function make_x_bounded( self )
      Spline2DMexWrapper( 'make_x_bounded', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Make spline surface computable outside the `x`-range is defined.
    %> If `x` is outside range value is extrapolated.
    %>
    %> ```{matlab}
    %>   spl.make_x_unbounded();
    %> ```
    %>
    function make_x_unbounded( self )
      Spline2DMexWrapper( 'make_x_unbounded', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Check if spline is computable only in the `x`-range where is defined.
    %> Return true if can be computed only in the range.
    %>
    %> ```{matlab}
    %>   ok = spl.is_x_bounded();
    %> ```
    %>
    function ok = is_x_bounded( self )
      ok = Spline2DMexWrapper( 'is_x_bounded', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Set spline surface as closed in `y direction
    %>
    %> ```{matlab}
    %>   spl.make_y_closed();
    %> ```
    %>
    function make_y_closed( self )
      Spline2DMexWrapper( 'make_y_closed', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Set spline surface as opened in `y direction
    %>
    %> ```{matlab}
    %>   spl.make_y_opened();
    %> ```
    %>
    function make_y_opened( self )
      Spline2DMexWrapper( 'make_y_opened', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Return true if spline surface is of closed type in `y`-direction
    %>
    %>
    %> ```{matlab}
    %>   ok = spl.is_y_closed();
    %> ```
    %>
    function ok = is_y_closed( self )
      ok = Spline2DMexWrapper( 'is_y_closed', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Make spline surface computable only in the `y`-range where is defined.
    %> If `y` is outside range an error is produced.
    %>
    %> ```{matlab}
    %>   spl.make_y_bounded();
    %> ```
    %>
    function make_y_bounded( self )
      Spline2DMexWrapper( 'make_y_bounded', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Make spline surface computable outside the `y`-range is defined.
    %> If `y` is outside range value is extrapolated.
    %>
    %> ```{matlab}
    %>   spl.make_y_unbounded();
    %> ```
    %>
    function make_y_unbounded( self )
      Spline2DMexWrapper( 'make_y_unbounded', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Check if spline is computable only in the `y`-range where is defined.
    %> Return true if can be computed only in the range.
    %>
    %> ```{matlab}
    %>   ok = spl.is_y_bounded();
    %> ```
    %>
    function ok = is_y_bounded( self )
      ok = Spline2DMexWrapper( 'is_y_bounded', self.objectHandle );
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
      res = Spline2DMexWrapper( 'x_min', self.objectHandle );
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
      res = Spline2DMexWrapper( 'x_max', self.objectHandle );
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
      res = Spline2DMexWrapper( 'y_min', self.objectHandle );
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
      res = Spline2DMexWrapper( 'y_max', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Get definition domain of the 2D-spline
    %>
    %> ```{matlab}
    %>   res = spl.z_min();
    %> ```
    %>
    function res = z_min( self )
      res = Spline2DMexWrapper( 'z_min', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Get definition domain of the 2D-spline
    %>
    %> ```{matlab}
    %>   res = spl.z_max();
    %> ```
    %>
    function res = z_max( self )
      res = Spline2DMexWrapper( 'z_max', self.objectHandle );
    end
  end
end
