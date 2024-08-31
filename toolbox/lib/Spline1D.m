  %>
  %> MATLAB class wrapper for the underlying C++ class
  %>
  classdef Spline1D < matlab.mixin.Copyable
  properties (SetAccess = private, Hidden = true)
    %> Handle to the underlying C++ class instance
    objectHandle;
    call_delete;
  end

  methods(Access = protected)
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %> Make a deep copy of a curve object
    %>
    %> **Usage**
    %>
    %> ```{matlab}
    %>   B = A.copy();
    %> ```
    %>
    %> where `A` is the curve object to be copied.
    %>
    function obj = copyElement( self )
      obj              = copyElement@matlab.mixin.Copyable(self);
      obj.objectHandle = feval( self.mexName, 'copy', self.objectHandle );
      obj.call_delete  = true;
    end
  end

  methods
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Build a spline given a table of point and type:
    %>
    %>
    %> ```{matlab}
    %>   % build a spline given points and spline spline
    %>   obj = Spline1D( x, y, kind );
    %>   % build Hermite type spline (do not use kind)
    %>   obj = Spline1D( x, y, yp );
    %> ```
    %>
    %> Kind is a string and can be any of
    %>
    %> |     kind     |   meaning                                                |
    %> | :----------- | :------------------------------------------------------- |
    %> | 'linear'     | linear spline (only continuous)                          |
    %> | 'cubic'      | cubic spline (\f$ C^2 \f$ function)                      |
    %> | 'akima'      | Akima non oscillatory spline (\f$ C^1 \f$ function)      |
    %> | 'bessel'     | bessel non oscillatory spline (\f$ C^1 \f$ function)     |
    %> | 'pchip'      | Monotone \f$ C^1 \f$ function                            |
    %> | 'hermite'    | Hermite spline (set \f$ p_k \f$ and \f$ p'_k \f$)        |
    %> | 'quintic'    | Quintic spline (\f$ C^3 \f$ function)                    |
    %>
    function self = Spline1D( kind, varargin )
      self.objectHandle = Spline1DMexWrapper( 'new', kind );
      obj.call_delete   = true;
      if nargin > 1
        Spline1DMexWrapper( 'build', self.objectHandle, varargin{:} );
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function delete(self)
      if self.call_delete
        Spline1DMexWrapper( 'delete', self.objectHandle );
      end
      self.objectHandle = 0;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Build a spline given a table of point and type:
    %>
    %>
    %> ```{matlab}
    %>   obj.build( x, y, kind );
    %>
    %>   obj.build( x, y, yp ); % build Hermite type spline
    %> ```
    %>
    function build( self, varargin )
      %> x, y [, yp or subtype]
      Spline1DMexWrapper( 'build', self.objectHandle, varargin{:} );
    end

    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Return the coeffs of the polynomial of the spline
    %>
    %>
    %> ```{matlab}
    %>   [coeffs,nodes] = obj.eval_coeffs();
    %> ```
    %>
    %> - [out] coeffs matrix which rows are the polynomial coeffs
    %> - [out] nodes where the spline is defined
    %>
    function varargout = eval_coeffs( self )
      [varargout{1:nargout}] = Spline1DMexWrapper( 'coeffs', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Search the max and min values of `y` along the spline
    %> with the corresponding `x` position
    %>
    %>
    %> ```{matlab}
    %>  [i_min_pos,x_min_pos,y_min,...
    %>   i_max_pos,x_max_pos,y_max] = obj.eval_y_min_max();
    %>
    %>  SS = obj.eval_y_min_max();
    %> ```
    %>
    %> - [out] i_min_pos where is the minimum (interval)
    %> - [out] x_min_pos where is the minimum
    %> - [out] y_min     the minimum value
    %> - [out] i_max_pos where is the maximum (interval)
    %> - [out] x_max_pos where is the maximum
    %> - [out] y_max     the maximum value
    %>
    %> - [out] SS.i_min_pos vector where is the local minimum (interval)
    %> - [out] SS.x_min_pos vector where is the local minimum
    %> - [out] SS.y_min     vector where is the local minimum value
    %> - [out] SS.i_max_pos vector where is the local maximum (interval)
    %> - [out] SS.x_max_pos vector where is the local maximum
    %> - [out] SS.y_max     vector where is the local maximum value
    %>
    function varargout = eval_y_min_max( self )
      [varargout{1:nargout}] = Spline1DMexWrapper( 'y_min_max', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Evaluate spline at `x` returning vaues and derivatives
    %>
    %> ```{matlab}
    %>   p = obj.eval( x );
    %>   [p,p_D] = obj.eval( x );
    %>   [p,p_D,p_DD] = obj.eval( x );
    %>   [p,p_D,p_DD,p_DDD] = obj.eval( x );
    %>   [p,p_D,p_DD,p_DDD,p_DDDD] = obj.eval( x );
    %>   [p,p_D],p_DD,p_DDD,p_DDDD,p_DDDDD] = obj.eval( x );
    %> ```
    %>
    function varargout = eval( self, x )
      varargout{1} = Spline1DMexWrapper( 'eval', self.objectHandle, x );
      if nargout >= 2
        varargout{2} = Spline1DMexWrapper( 'eval_D', self.objectHandle, x );
      end
      if nargout >= 3
        varargout{3} = Spline1DMexWrapper( 'eval_DD', self.objectHandle, x );
      end
      if nargout >= 4
        varargout{4} = Spline1DMexWrapper( 'eval_DDD', self.objectHandle, x );
      end
      if nargout >= 5
        varargout{5} = Spline1DMexWrapper( 'eval_DDDD', self.objectHandle, x );
      end
      if ~( nargout == 1 || nargout == 2 || nargout == 3 )
        error( 'Spline1D.eval, nargout = %d must be 1, 2, 3, 4 or 5\n', nargout );
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Evaluate spline derivative at `x`
    %>
    %>
    %> ```{matlab}
    %>   p_D = obj.eval_D( x );
    %> ```
    %>
    function dp = eval_D( self, x )
      dp = Spline1DMexWrapper( 'eval_D', self.objectHandle, x );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Evaluate spline second derivative at `x`
    %>
    %> ```{matlab}
    %>   p_DD = obj.eval_DD( x );
    %> ```
    %>
    function ddp = eval_DD( self, x )
      ddp = Spline1DMexWrapper( 'eval_DD', self.objectHandle, x );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Evaluate spline third derivative at `x`
    %>
    %> ```{matlab}
    %>   p_DDD = obj.eval_DDD( x );
    %> ```
    %>
    function dddp = eval_DDD( self, x )
      dddp = Spline1DMexWrapper( 'eval_DDD', self.objectHandle, x );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Evaluate spline 4th derivative at `x`
    %>
    %> ```{matlab}
    %>   p_DDDD = obj.eval_DDDD( x );
    %> ```
    %>
    function ddddp = eval_DDDD( self, x )
      ddddp = Spline1DMexWrapper( 'eval_DDDD', self.objectHandle, x );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Evaluate spline 5th derivative at `x`
    %>
    %> ```{matlab}
    %>   p_DDDDD = obj.eval_DDDDD( x );
    %> ```
    %>
    function dddddp = eval_DDDDD( self, x )
      dddddp = Spline1DMexWrapper( 'eval_DDDDD', self.objectHandle, x );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Set spline as closed
    %>
    %> ```{matlab}
    %>   obj.make_closed();
    %> ```
    %>
    function make_closed( self )
      Spline1DMexWrapper( 'make_closed', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Set spline as opened
    %>
    %> ```{matlab}
    %>   obj.make_opened();
    %> ```
    %>
    function make_opened( self )
      Spline1DMexWrapper( 'make_opened', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Return true if spline is of closed type
    %>
    %> ```{matlab}
    %>   ok = obj.is_closed();
    %> ```
    %>
    function ok = is_closed( self )
      ok = Spline1DMexWrapper( 'is_closed', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Make spline computable only in the range is defined.
    %> If `x` is outside range an error is produced.
    %>
    %> ```{matlab}
    %>   obj.make_bounded();
    %> ```
    %>
    function make_bounded( self )
      Spline1DMexWrapper( 'make_bounded', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Make spline computable outside the range is defined.
    %> If `x` is outside range value is extrapolated.
    %>
    %> ```{matlab}
    %>   obj.make_unbounded();
    %> ```
    %>
    function make_unbounded( self )
      Spline1DMexWrapper( 'make_unbounded', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Check if spline is computable only in the range is defined.
    %> Return true if can be computed only in the range.
    %>
    %> ```{matlab}
    %>   ok = obj.is_bounded();
    %> ```
    %>
    function ok = is_bounded( self )
      ok = Spline1DMexWrapper( 'is_bounded', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Set extrapolation outside the range as a constant spline
    %>
    %> ```{matlab}
    %>   obj.make_extended_constant();
    %> ```
    %>
    function make_extended_constant( self )
      Spline1DMexWrapper( 'make_extended_constant', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Set extrapolation outside the range as a polynomial
    %>
    %> ```{matlab}
    %>   obj.make_extended_not_constant();
    %> ```
    %>
    function make_extended_not_constant( self )
      Spline1DMexWrapper( 'make_extended_not_constant', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Check if extrapolation outside the range as a polynomial
    %>
    %> ```{matlab}
    %>   ok = obj.is_extended_constant();
    %> ```
    %>
    function ok = is_extended_constant( self )
      ok = Spline1DMexWrapper( 'is_extended_constant', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Return initial x-coordinate of the spline
    %>
    %> ```{matlab}
    %>   x = obj.x_begin();
    %> ```
    %>
    function x = x_begin( self )
      x = Spline1DMexWrapper( 'x_begin', self.objectHandle );
    end
    %>
    %> \deprecated will be removed in future version
    %>
    function x = xBegin( self )
      x = Spline1DMexWrapper( 'x_begin', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Return final x-coordinate of the spline
    %>
    %> ```{matlab}
    %>   x = obj.x_end();
    %> ```
    %>
    function x = xEnd( self )
      x = Spline1DMexWrapper( 'x_end', self.objectHandle );
    end
    %>
    %> \deprecated will be removed in future version
    %>
    function x = x_end( self )
      x = Spline1DMexWrapper( 'x_end', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Return initial y-coordinate of the spline
    %>
    %>
    %> ```{matlab}
    %>   y = obj.y_begin();
    %> ```
    %>
    function y = y_begin( self )
      y = Spline1DMexWrapper( 'y_begin', self.objectHandle );
    end
    %>
    %> \deprecated will be removed in future version
    %>
    function y = yBegin( self )
      y = Spline1DMexWrapper( 'y_begin', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Return final y-coordinate of the spline
    %>
    %> ```{matlab}
    %>   y = obj.y_end();
    %> ```
    %>
    function y = y_end( self )
      y = Spline1DMexWrapper( 'y_end', self.objectHandle );
    end
    %>
    %> \deprecated will be removed in future version
    %>
    function y = yEnd( self )
      y = Spline1DMexWrapper( 'y_end', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Return minimum x-coordinate of the spline
    %>
    %> ```{matlab}
    %>   x = obj.x_min();
    %> ```
    %>
    function x = x_min( self )
      x = Spline1DMexWrapper( 'x_min', self.objectHandle );
    end
    %>
    %> \deprecated will be removed in future version
    %>
    function x = xMin( self )
      x = Spline1DMexWrapper( 'x_min', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Return maximum x-coordinate of the spline
    %>
    %> ```{matlab}
    %>   x = obj.x_max();
    %> ```
    %>
    function x = x_max( self )
      x = Spline1DMexWrapper( 'x_max', self.objectHandle );
    end
    %>
    %> \deprecated will be removed in future version
    %>
    function x = xMax( self )
      x = Spline1DMexWrapper( 'x_max', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Return minimum y-coordinate of the spline
    %>
    %> ```{matlab}
    %>   y = obj.y_min();
    %> ```
    %>
    function y = y_min( self )
      y = Spline1DMexWrapper( 'y_min', self.objectHandle );
    end
    %>
    %> \deprecated will be removed in future version
    %>
    function y = yMin( self )
      y = Spline1DMexWrapper( 'y_min', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Return maximum y-coordinate of the spline
    %>
    %> ```{matlab}
    %>   y = obj.y_max();
    %> ```
    %>
    function y = y_max( self )
      y = Spline1DMexWrapper( 'y_max', self.objectHandle );
    end
    %>
    %> \deprecated will be removed in future version
    %>
    function y = yMax( self )
      y = Spline1DMexWrapper( 'y_max', self.objectHandle );
    end
  end
end
