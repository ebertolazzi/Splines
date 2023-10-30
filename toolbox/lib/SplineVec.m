%>
%> MATLAB class wrapper for the underlying C++ class
%>
%> The construction of the spline is done as follows
%>
%> - instantiate the spline object
%>
%> \rst
%>
%>   .. code-block:: matlab
%>
%>     spl = SplineVec();
%>
%> \endrst
%>
%> - load the points as a matrix `dim x npts` where
%>   `dim` is the space dimension and `npts` is the
%>    number of poinst
%>
%> \rst
%>
%>   .. code-block:: matlab
%>
%>     spl.setup( P ); % P is a dim x npts matrix
%>
%> \endrst
%>
%> - build or load the knots
%>
%> \rst
%>
%>   .. code-block:: matlab
%>
%>     spl.knots( T );   % vector of knots
%>     spl.chordal();    % build knots using chordal distance
%>     spl.centripetal(; % build knots centripetal
%>
%> \endrst
%>
%> - build the spline, for the moment only `Catmull Rom`
%>
%> \rst
%>
%>   .. code-block:: matlab
%>
%>     spl.CatmullRom();
%>
%> \endrst
%>
classdef SplineVec < matlab.mixin.Copyable
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
    %> \rst
    %> .. code-block:: matlab
    %>
    %>   B = A.copy();
    %>
    %> \endrst
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
    %> Build an empty parametric spline:
    %>
    %> \rst
    %>
    %>   .. code-block:: matlab
    %>
    %>     spl = SplineVec();
    %>
    %> \endrst
    %>
    function self = SplineVec()
      self.objectHandle = SplineVecMexWrapper( 'new' );
      obj.call_delete   = true;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function delete(self)
      if self.call_delete
        SplineVecMexWrapper( 'delete', self.objectHandle );
      end
      self.objectHandle = 0;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Setup the point of the spline
    %>
    %> \rst
    %>
    %>   .. code-block:: matlab
    %>
    %>     spl.setup( P );
    %>
    %> \endrst
    %>
    function setup( self, P )
      SplineVecMexWrapper( 'setup', self.objectHandle, P );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Setup the knots of the spline
    %>
    %> \rst
    %>
    %>   .. code-block:: matlab
    %>
    %>     spl.knots( t );
    %>
    %> \endrst
    %>
    function knots( self, t )
      SplineVecMexWrapper( 'knots', self.objectHandle, t );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Get the knots of the spline
    %>
    %> \rst
    %>
    %>   .. code-block:: matlab
    %>
    %>     t = spl.get_knots();
    %>
    %> \endrst
    %>
    function k = get_knots( self )
      k = SplineVecMexWrapper( 'get_knots', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Build the knots using chordal distance
    %>
    %> \rst
    %>
    %>   .. math::
    %>
    %>     t_0 = 0, \qquad t_k = t_{k-1} + || \mathbf{p}_k - \mathbf{p}_{k-1} ||
    %>
    %> \endrst
    %>
    function chordal( self )
      SplineVecMexWrapper( 'chordal', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Build the knots using centripetal parametrization
    %>
    %> \rst
    %>
    %>   .. math::
    %>
    %>     t_0 = 0, \qquad t_k = t_{k-1} + || \mathbf{p}_k - \mathbf{p}_{k-1} ||^{1/2}
    %>
    %> \endrst
    %>
    function centripetal( self )
      SplineVecMexWrapper( 'centripetal', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Build the spline using Catmull Rom algorithm
    %>
    %> \rst
    %>
    %>   .. code-block:: matlab
    %>
    %>     P = [ 0, 0; 1.34, 5; 5, 8.66; 10, 10; 10.6, 10.4; 10.7, 12; ...
    %>           10.7, 28.6; 10.8, 30.2; 11.4, 30.6; 19.6, 30.6; ...
    %>           20.2, 30.2; 20.3, 28.6; 20.3, 12; 20.4, 10.4; ...
    %>           21, 10; 26, 8.66; 29.66, 5; 31, 0 ];
    %>
    %>     S = SplineVec();  % initialize spline
    %>     S.setup(P.');     % load points
    %>     S.centripetal();  % setup knots
    %>
    %>     % plot interpolation points
    %>     hold off;
    %>     plot( P(:,1), P(:,2), 'o','Color','red', ...
    %>           'MarkerSize',10,'MarkerFaceColor','blue','MarkerEdgeColor','green');
    %>     hold on;
    %>     % sample spline for plotting
    %>     t  = linspace(S.tmin(),S.tmax(),1000);
    %>     PP = S.eval(t);
    %>
    %>     % plot curve
    %>     plot( PP(1,:), PP(2,:), '-', 'Color', 'blue', 'Linewidth', 3);
    %>
    %>   .. image:: ../../images/exampleVec.png
    %>      :width: 80%
    %>      :align: center
    %>
    %> \endrst
    %>
    function CatmullRom( self )
      SplineVecMexWrapper( 'CatmullRom', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Evaluate spline at `x`
    %>
    %> \rst
    %>
    %>   .. code-block:: matlab
    %>
    %>     p = obj.eval( x );
    %>
    %> \endrst
    %>
    function p = eval( self, x )
      p = SplineVecMexWrapper( 'eval', self.objectHandle, x );
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
      dp = SplineVecMexWrapper( 'eval_D', self.objectHandle, x );
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
      ddp = SplineVecMexWrapper( 'eval_DD', self.objectHandle, x );
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
      dddp = SplineVecMexWrapper( 'eval_DDD', self.objectHandle, x );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Evaluate spline curvature at `x`
    %>
    %> \rst
    %>
    %>   .. code-block:: matlab
    %>
    %>     kappa = obj.curvature( x );
    %>
    %> \endrst
    %>
    function cur = curvature( self, x )
      cur = SplineVecMexWrapper( 'eval_curvature', self.objectHandle, x );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Evaluate spline curvature derivative at `x`
    %>
    %> \rst
    %>
    %>   .. code-block:: matlab
    %>
    %>     kappa = obj.curvature_D( x );
    %>
    %> \endrst
    %>
    function cur_D = curvature_D( self, x )
      cur_D = SplineVecMexWrapper( 'eval_curvature_D', self.objectHandle, x );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Return initial t-coordinate of the spline
    %>
    %> \rst
    %>
    %>   .. code-block:: matlab
    %>
    %>     t = obj.tmin();
    %>
    %> \endrst
    %>
    function t = tmin( self, x )
      t = SplineVecMexWrapper( 'tmin', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %>
    %> Return initial t-coordinate of the spline
    %>
    %> \rst
    %>
    %>   .. code-block:: matlab
    %>
    %>     t = obj.tmax();
    %>
    %> \endrst
    %>
    function t = tmax( self, x )
      t = SplineVecMexWrapper( 'tmax', self.objectHandle );
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end
end
