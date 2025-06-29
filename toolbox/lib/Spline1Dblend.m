  %>
  %> MATLAB class wrapper for the underlying C++ class
  %>
  classdef Spline1Dblend < handle
  
  properties (SetAccess = private, Hidden = true)
    %> Handle to the underlying C++ class instance
    objectHandle;
    call_delete;
  end

  methods
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Build a spline given a table of point and type:
    %>
    %>
    %> ```{matlab}
    %>   % build a spline given points and spline spline
    %>   obj = Spline1Dblend( x, y, kind );
    %>   % build Hermite type spline (do not use kind)
    %>   obj = Spline1Dblend( x, y, yp );
    %> ```
    %>
    %> Kind is a string and can be any of
    %>
    %> |     kind         |   meaning                                                |
    %> | :-----------     | :------------------------------------------------------- |
    %> | 'linear'         | linear spline (only continuous)                          |
    %> | 'cubic'          | cubic spline (\f$ C^2 \f$ function)                      |
    %> | 'akima'          | Akima non oscillatory spline (\f$ C^1 \f$ function)      |
    %> | 'bessel'         | bessel non oscillatory spline (\f$ C^1 \f$ function)     |
    %> | 'pchip'          | Monotone \f$ C^1 \f$ function                            |
    %> | 'hermite'        | Hermite spline (set \f$ p_k \f$ and \f$ p'_k \f$)        |
    %> | 'quintic'        | Quintic spline (\f$ C^3 \f$ function)                    |
    %> | 'quintic_pchip'  | Quintic spline (\f$ C^2 \f$ function)                    |
    %> | 'quintic_akima'  | Quintic spline (\f$ C^2 \f$ function)                    |
    %> | 'quintic_bessel' | Quintic spline (\f$ C^2 \f$ function)                    |
    %>
    function self = Spline1Dblend( S )
      self.objectHandle = Spline1DblendMexWrapper( 'new', S );
      self.call_delete  = true;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function delete(self)
      if self.call_delete
        Spline1DblendMexWrapper( 'delete', self.objectHandle );
      end
      self.call_delete  = false;
      self.objectHandle = 0;
    end

    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Evaluate spline at `x` parameter `s` returning vaues and derivatives
    %>
    %> ```{matlab}
    %>   p = spl.eval( x, s );
    %>   [p,p_D] = spl.eval( x, s );
    %>   [p,p_D,p_DD] = spl.eval( x, s );
    %>   [p,p_D,p_DD,p_DDD] = spl.eval( x, s );
    %>   [p,p_D,p_DD,p_DDD,p_DDDD] = spl.eval( x, s );
    %>   [p,p_D],p_DD,p_DDD,p_DDDD,p_DDDDD] = spl.eval( x, s );
    %> ```
    %>
    function varargout = eval( self, x, s )
      varargout{1} = Spline1DblendMexWrapper( 'eval', self.objectHandle, x, s );
      if nargout >= 2
        varargout{2} = Spline1DblendMexWrapper( 'eval_D', self.objectHandle, x, s );
      end
      if nargout >= 3
        varargout{3} = Spline1DblendMexWrapper( 'eval_DD', self.objectHandle, x, s );
      end
      if nargout >= 4
        varargout{4} = Spline1DblendMexWrapper( 'eval_DDD', self.objectHandle, x, s );
      end
      if nargout >= 5
        varargout{5} = Spline1DblendMexWrapper( 'eval_DDDD', self.objectHandle, x, s );
      end
      if ~( nargout == 1 || nargout == 2 || nargout == 3 )
        error( 'Spline1Dblend.eval, nargout = %d must be 1, 2, 3, 4 or 5\n', nargout );
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Evaluate spline derivative at `x` parameter `s`
    %>
    %>
    %> ```{matlab}
    %>   p_D = spl.eval_D( x, s );
    %> ```
    %>
    function dp = eval_D( self, x )
      dp = Spline1DblendMexWrapper( 'eval_D', self.objectHandle, x, s );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Evaluate spline second derivative at `x` parameter `s`
    %>
    %> ```{matlab}
    %>   p_DD = spl.eval_DD( x, s );
    %> ```
    %>
    function ddp = eval_DD( self, x, s )
      ddp = Spline1DblendMexWrapper( 'eval_DD', self.objectHandle, x, s );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Evaluate spline third derivative at `x` parameter `s`
    %>
    %> ```{matlab}
    %>   p_DDD = spl.eval_DDD( x, s );
    %> ```
    %>
    function dddp = eval_DDD( self, x, s )
      dddp = Spline1DblendMexWrapper( 'eval_DDD', self.objectHandle, x, s );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Evaluate spline 4th derivative at `x` parameter `s`
    %>
    %> ```{matlab}
    %>   p_DDDD = spl.eval_DDDD( x, s );
    %> ```
    %>
    function ddddp = eval_DDDD( self, x, s )
      ddddp = Spline1DblendMexWrapper( 'eval_DDDD', self.objectHandle, x, s );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Evaluate spline 5th derivative at `x` parameter `s`
    %>
    %> ```{matlab}
    %>   p_DDDDD = spl.eval_DDDDD( x, s );
    %> ```
    %>
    function dddddp = eval_DDDDD( self, x )
      dddddp = Spline1DblendMexWrapper( 'eval_DDDDD', self.objectHandle, x, s );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Return initial x-coordinate of the spline
    %>
    %> ```{matlab}
    %>   x = spl.x_begin( s );
    %> ```
    %>
    function x = x_begin( self, s )
      x = Spline1DblendMexWrapper( 'x_begin', self.objectHandle, s );
    end
    %>
    %> \deprecated will be removed in future version
    %>
    function x = xBegin( self, s )
      x = Spline1DblendMexWrapper( 'x_begin', self.objectHandle, s );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Return final x-coordinate of the spline
    %>
    %> ```{matlab}
    %>   x = spl.x_end( s );
    %> ```
    %>
    function x = x_end( self, s )
      x = Spline1DblendMexWrapper( 'x_end', self.objectHandle, s );
    end
    %>
    %> \deprecated will be removed in future version
    %>
    function x = xEnd( self, s )
      x = Spline1DblendMexWrapper( 'x_end', self.objectHandle, s );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Return initial y-coordinate of the spline
    %>
    %>
    %> ```{matlab}
    %>   y = spl.y_begin( s );
    %> ```
    %>
    function y = y_begin( self, s )
      y = Spline1DblendMexWrapper( 'y_begin', self.objectHandle, s );
    end
    %>
    %> \deprecated will be removed in future version
    %>
    function y = yBegin( self )
      y = Spline1DblendMexWrapper( 'y_begin', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Return final y-coordinate of the spline
    %>
    %> ```{matlab}
    %>   y = spl.y_end( s );
    %> ```
    %>
    function y = y_end( self, s )
      y = Spline1DblendMexWrapper( 'y_end', self.objectHandle, s );
    end
    %>
    %> \deprecated will be removed in future version
    %>
    function y = yEnd( self, s )
      y = Spline1DblendMexWrapper( 'y_end', self.objectHandle, s );
    end
  end
end
