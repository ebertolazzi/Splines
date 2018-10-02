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
    function p = eval( self, x, y )
      p = Spline2DMexWrapper( 'eval', self.objectHandle, x, y );
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
  end
end
