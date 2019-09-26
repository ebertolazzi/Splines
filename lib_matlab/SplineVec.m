classdef SplineVec < handle
  %% MATLAB class wrapper for the underlying C++ class
  properties (SetAccess = private, Hidden = true)
    objectHandle; % Handle to the underlying C++ class instance
  end

  methods
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = SplineVec()
      self.objectHandle = SplineVecMexWrapper( 'new' );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function delete(self)
      % Destroy the C++ class instance
      SplineVecMexWrapper( 'delete', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function setup( self, y )
      SplineVecMexWrapper( 'setup', self.objectHandle, y );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function knots( self, x )
      SplineVecMexWrapper( 'knots', self.objectHandle, x );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function k = get_knots( self )
      k = SplineVecMexWrapper( 'getNodes', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function chord( self )
      SplineVecMexWrapper( 'chord', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function centripetal( self )
      SplineVecMexWrapper( 'centripetal', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function CatmullRom( self )
      SplineVecMexWrapper( 'CatmullRom', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function p = eval( self, x )
      p = SplineVecMexWrapper( 'eval', self.objectHandle, x );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function dp = eval_D( self, x )
      dp = SplineVecMexWrapper( 'eval_D', self.objectHandle, x );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function ddp = eval_DD( self, x )
      ddp = SplineVecMexWrapper( 'eval_DD', self.objectHandle, x );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function dddp = eval_DDD( self, x )
      dddp = SplineVecMexWrapper( 'eval_DDD', self.objectHandle, x );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function cur = curvature( self, x )
      cur = SplineVecMexWrapper( 'eval_curvature', self.objectHandle, x );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function cur_D = curvature_D( self, x )
      cur_D = SplineVecMexWrapper( 'eval_curvature_D', self.objectHandle, x );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function t = tmin( self, x )
      t = SplineVecMexWrapper( 'tmin', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function t = tmax( self, x )
      t = SplineVecMexWrapper( 'tmax', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end
end
