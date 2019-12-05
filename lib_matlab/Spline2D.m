classdef Spline2D < handle
  %% MATLAB class wrapper for the underlying C++ class
  properties (SetAccess = private, Hidden = true)
    objectHandle; % Handle to the underlying C++ class instance
  end

  methods
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = Spline2D( name, varargin )
      self.objectHandle = Spline2DMexWrapper( 'new', name );
      if nargin > 1
        Spline2DMexWrapper( 'build', self.objectHandle, ...
                             varargin{1}, varargin{2}, varargin{3} );
      end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function delete(self)
      % Destroy the C++ class instance
      Spline2DMexWrapper( 'delete', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function build( self, x, y, z )
      Spline2DMexWrapper( 'build', self.objectHandle, x, y );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
    function Dx = eval_Dx( self, x, y )
      Dx = Spline2DMexWrapper( 'eval_Dx', self.objectHandle, x, y );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function Dy = eval_Dy( self, x, y )
      Dy = Spline2DMexWrapper( 'eval_Dy', self.objectHandle, x, y );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function Dxx = eval_Dxx( self, x, y )
      Dxx = Spline2DMexWrapper( 'eval_Dxx', self.objectHandle, x, y );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function Dxy = eval_Dxy( self, x, y )
      Dxy = Spline2DMexWrapper( 'eval_Dxy', self.objectHandle, x, y );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function Dyy = eval_Dyy( self, x, y )
      Dyy = Spline2DMexWrapper( 'eval_Dyy', self.objectHandle, x, y );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function make_x_closed( self )
      Spline2DMexWrapper( 'make_x_closed', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function make_x_opened( self )
      Spline2DMexWrapper( 'make_x_opened', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function ok = is_x_closed( self )
      ok = Spline2DMexWrapper( 'is_x_closed', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function make_x_bounded( self )
      Spline2DMexWrapper( 'make_x_bounded', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function make_x_unbounded( self )
      Spline2DMexWrapper( 'make_x_unbounded', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function ok = is_x_bounded( self )
      ok = Spline2DMexWrapper( 'is_x_bounded', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function make_y_closed( self )
      Spline2DMexWrapper( 'make_y_closed', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function make_y_opened( self )
      Spline2DMexWrapper( 'make_y_opened', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function ok = is_y_closed( self )
      ok = Spline2DMexWrapper( 'is_y_closed', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function make_y_bounded( self )
      Spline2DMexWrapper( 'make_y_bounded', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function make_y_unbounded( self )
      Spline2DMexWrapper( 'make_y_unbounded', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function ok = is_y_bounded( self )
      ok = Spline2DMexWrapper( 'is_y_bounded', self.objectHandle );
    end
  end
end
