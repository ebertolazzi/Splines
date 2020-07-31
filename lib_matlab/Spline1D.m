classdef Spline1D < handle
  %% MATLAB class wrapper for the underlying C++ class
  properties (SetAccess = private, Hidden = true)
    objectHandle; % Handle to the underlying C++ class instance
  end

  methods
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = Spline1D( kind, varargin )
      % kind [, t, pnts, subtype ]
      self.objectHandle = Spline1DMexWrapper( 'new', kind );
      if nargin > 1
        Spline1DMexWrapper( 'build', self.objectHandle, varargin{:} );
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function delete(self)
      % Destroy the C++ class instance
      Spline1DMexWrapper( 'delete', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function build( self, varargin )
      % x, y [, yp or subtype]
      Spline1DMexWrapper( 'build', self.objectHandle, varargin{:} );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
    function dp = eval_D( self, x )
      dp = Spline1DMexWrapper( 'eval_D', self.objectHandle, x );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function ddp = eval_DD( self, x )
      ddp = Spline1DMexWrapper( 'eval_DD', self.objectHandle, x );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function dddp = eval_DDD( self, x )
      dddp = Spline1DMexWrapper( 'eval_DDD', self.objectHandle, x );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function ddddp = eval_DDDD( self, x )
      ddddp = Spline1DMexWrapper( 'eval_DDDD', self.objectHandle, x );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function dddddp = eval_DDDDD( self, x )
      dddddp = Spline1DMexWrapper( 'eval_DDDDD', self.objectHandle, x );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function make_closed( self )
      Spline1DMexWrapper( 'make_closed', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function make_opened( self )
      Spline1DMexWrapper( 'make_opened', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function ok = is_closed( self )
      ok = Spline1DMexWrapper( 'is_closed', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function make_bounded( self )
      Spline1DMexWrapper( 'make_bounded', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function make_unbounded( self )
      Spline1DMexWrapper( 'make_unbounded', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function ok = is_bounded( self )
      ok = Spline1DMexWrapper( 'is_bounded', self.objectHandle );
    end
  end
end
