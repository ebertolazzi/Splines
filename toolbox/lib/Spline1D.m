  %> MATLAB class wrapper for the underlying C++ class
  classdef Spline1D < handle
  properties (SetAccess = private, Hidden = true)
    %> Handle to the underlying C++ class instance
    objectHandle;
  end

  methods
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Build a spline given a table of point and type:
    %>
    %> \rst
    %>
    %>   .. code-block:: matlab
    %>
    %>     obj = Spline1D( x, y, kind );
    %>
    %>     obj = Spline1D( x, y, yp ); % build Hermite type spline, not use kind
    %>
    %> \endrst
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
      %> kind [, t, pnts, subtype ]
      self.objectHandle = Spline1DMexWrapper( 'new', kind );
      if nargin > 1
        Spline1DMexWrapper( 'build', self.objectHandle, varargin{:} );
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function delete(self)
      Spline1DMexWrapper( 'delete', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Build a spline given a table of point and type:
    %>
    %> \rst
    %>
    %>   .. code-block:: matlab
    %>
    %>     obj.build( x, y, kind );
    %>
    %>     obj.build( x, y, yp ); % build Hermite type spline
    %>
    %> \endrst
    %>
    function build( self, varargin )
      %> x, y [, yp or subtype]
      Spline1DMexWrapper( 'build', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Evaluate spline at `x` returning vaues and derivatives
    %>
    %> \rst
    %>
    %>   .. code-block:: matlab
    %>
    %>     p = obj.eval( x );
    %>     [p,p_D] = obj.eval( x );
    %>     [p,p_D,p_DD] = obj.eval( x );
    %>     [p,p_D,p_DD,p_DDD] = obj.eval( x );
    %>     [p,p_D,p_DD,p_DDD,p_DDDD] = obj.eval( x );
    %>     [p,p_D],p_DD,p_DDD,p_DDDD,p_DDDDD] = obj.eval( x );
    %>
    %> \endrst
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
    %> \rst
    %>
    %>   .. code-block:: matlab
    %>
    %>     p_D = obj.eval_D( x );
    %>
    %> \endrst
    %>
    function dp = eval_D( self, x )
      dp = Spline1DMexWrapper( 'eval_D', self.objectHandle, x );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Evaluate spline second derivative at `x`
    %>
    %> \rst
    %>
    %>   .. code-block:: matlab
    %>
    %>     p_DD = obj.eval_DD( x );
    %>
    %> \endrst
    %>
    function ddp = eval_DD( self, x )
      ddp = Spline1DMexWrapper( 'eval_DD', self.objectHandle, x );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Evaluate spline third derivative at `x`
    %>
    %> \rst
    %>
    %>   .. code-block:: matlab
    %>
    %>     p_DDD = obj.eval_DDD( x );
    %>
    %> \endrst
    %>
    function dddp = eval_DDD( self, x )
      dddp = Spline1DMexWrapper( 'eval_DDD', self.objectHandle, x );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Evaluate spline 4th derivative at `x`
    %>
    %> \rst
    %>
    %>   .. code-block:: matlab
    %>
    %>     p_DDDD = obj.eval_DDDD( x );
    %>
    %> \endrst
    %>
    function ddddp = eval_DDDD( self, x )
      ddddp = Spline1DMexWrapper( 'eval_DDDD', self.objectHandle, x );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Evaluate spline 5th derivative at `x`
    %>
    %> \rst
    %>
    %>   .. code-block:: matlab
    %>
    %>     p_DDDDD = obj.eval_DDDDD( x );
    %>
    %> \endrst
    %>
    function dddddp = eval_DDDDD( self, x )
      dddddp = Spline1DMexWrapper( 'eval_DDDDD', self.objectHandle, x );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Set spline as closed
    %>
    %> \rst
    %>
    %>   .. code-block:: matlab
    %>
    %>     obj.make_closed();
    %>
    %> \endrst
    %>
    function make_closed( self )
      Spline1DMexWrapper( 'make_closed', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Set spline as opened
    %>
    %> \rst
    %>
    %>   .. code-block:: matlab
    %>
    %>     obj.make_opened();
    %>
    %> \endrst
    %>
    function make_opened( self )
      Spline1DMexWrapper( 'make_opened', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Return true if spline is of closed type
    %>
    %> \rst
    %>
    %>   .. code-block:: matlab
    %>
    %>     ok = obj.is_closed();
    %>
    %> \endrst
    %>
    function ok = is_closed( self )
      ok = Spline1DMexWrapper( 'is_closed', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Make spline computable only in the range is defined.
    %> If `x` is outside range an error is produced.
    %>
    %> \rst
    %>
    %>   .. code-block:: matlab
    %>
    %>     obj.make_bounded();
    %>
    %> \endrst
    %>
    function make_bounded( self )
      Spline1DMexWrapper( 'make_bounded', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Make spline computable outside the range is defined.
    %> If `x` is outside range value is extrapolated.
    %>
    %> \rst
    %>
    %>   .. code-block:: matlab
    %>
    %>     obj.make_unbounded();
    %>
    %> \endrst
    %>
    function make_unbounded( self )
      Spline1DMexWrapper( 'make_unbounded', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Check if spline is computable only in the range is defined.
    %> Return true if can be computed only in the range.
    %>
    %> \rst
    %>
    %>   .. code-block:: matlab
    %>
    %>     ok = obj.is_bounded();
    %>
    %> \endrst
    %>
    function ok = is_bounded( self )
      ok = Spline1DMexWrapper( 'is_bounded', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Set extrapolation outside the range as a constant spline
    %>
    %> \rst
    %>
    %>   .. code-block:: matlab
    %>
    %>     obj.make_extended_constant();
    %>
    %> \endrst
    %>
    function make_extended_constant( self )
      Spline1DMexWrapper( 'make_extended_constant', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Set extrapolation outside the range as a polynomial
    %>
    %> \rst
    %>
    %>   .. code-block:: matlab
    %>
    %>     obj.make_extended_not_constant();
    %>
    %> \endrst
    %>
    function make_extended_not_constant( self )
      Spline1DMexWrapper( 'make_extended_not_constant', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Check if extrapolation outside the range as a polynomial
    %>
    %> \rst
    %>
    %>   .. code-block:: matlab
    %>
    %>     ok = obj.is_extended_constant();
    %>
    %> \endrst
    %>
    function ok = is_extended_constant( self )
      ok = Spline1DMexWrapper( 'is_extended_constant', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Return initial x-coordinate of the spline
    %>
    %> \rst
    %>
    %>   .. code-block:: matlab
    %>
    %>     x = obj.xBegin();
    %>
    %> \endrst
    %>
    function x = xBegin( self )
      x = Spline1DMexWrapper( 'xBegin', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Return final x-coordinate of the spline
    %>
    %> \rst
    %>
    %>   .. code-block:: matlab
    %>
    %>     x = obj.xEnd();
    %>
    %> \endrst
    %>
    function x = xEnd( self )
      x = Spline1DMexWrapper( 'xEnd', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Return initial y-coordinate of the spline
    %>
    %> \rst
    %>
    %>   .. code-block:: matlab
    %>
    %>     y = obj.yBegin();
    %>
    %> \endrst
    %>
    function y = yBegin( self )
      y = Spline1DMexWrapper( 'yBegin', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Return final y-coordinate of the spline
    %>
    %> \rst
    %>
    %>   .. code-block:: matlab
    %>
    %>     y = obj.yEnd();
    %>
    %> \endrst
    %>
    function y = yEnd( self )
      y = Spline1DMexWrapper( 'yEnd', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Return minimum x-coordinate of the spline
    %>
    %> \rst
    %>
    %>   .. code-block:: matlab
    %>
    %>     x = obj.xMin();
    %>
    %> \endrst
    %>
    function x = xMin( self )
      x = Spline1DMexWrapper( 'xMin', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Return maximum x-coordinate of the spline
    %>
    %> \rst
    %>
    %>   .. code-block:: matlab
    %>
    %>     x = obj.xMax();
    %>
    %> \endrst
    %>
    function x = xMax( self )
      x = Spline1DMexWrapper( 'xMax', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Return minimum y-coordinate of the spline
    %>
    %> \rst
    %>
    %>   .. code-block:: matlab
    %>
    %>     y = obj.yMin();
    %>
    %> \endrst
    %>
    function y = yMin( self )
      y = Spline1DMexWrapper( 'yMin', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Return maximum y-coordinate of the spline
    %>
    %> \rst
    %>
    %>   .. code-block:: matlab
    %>
    %>     y = obj.yMax();
    %>
    %> \endrst
    %>
    function y = yMax( self )
      y = Spline1DMexWrapper( 'yMax', self.objectHandle );
    end
  end
end
