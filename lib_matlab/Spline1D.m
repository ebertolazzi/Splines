classdef Spline1D < handle
  %% MATLAB class wrapper for the underlying C++ class
  properties (SetAccess = private, Hidden = true)
    objectHandle; % Handle to the underlying C++ class instance
  end

  methods
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = Spline1D( kind, varargin )
      self.objectHandle = Spline1DMexWrapper( 'new', kind );
      if nargin > 1
        Spline1DMexWrapper( 'build', self.objectHandle, varargin{1}, varargin{2} );
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function delete(self)
      % Destroy the C++ class instance
      Spline1DMexWrapper( 'delete', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function build( self, x, y )
      Spline1DMexWrapper( 'build', self.objectHandle, x, y );
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
      if ~( nargout == 1 || nargout == 2 || nargout == 3 ) 
        error( 'Spline2D.eval, nargout = %d must be 1, 2 or 3\n', nargout);
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
  end
end
