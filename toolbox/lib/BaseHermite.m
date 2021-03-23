classdef BaseHermite < handle

  methods

    function self = BaseHermite()
      %
      % Build a matlab object storing Hermite base evaluator
      %
    end

    function delete(self)
      %
      % Destroy a matlab object storing Hermite base evaluator
      %
    end
    %
    % --------------------------------------------------------------------
    %
    function varargout = base( ~, varargin )
      %
      %  Evaluate an Hermite base (cubic degree) at point(s) t for the interval 
      %  :math:`[0,1]` (or :math:`[0,H]` if second argument
      %  is present).
      %
      %  .. code-block:: matlab
      %
      %      BASE = object.base( t )      % base for the interval [0,1]
      %      BASE = object.base( t, H )   % base for the interval [0,H]
      %
      %  BASE is a matrix length(t) x 4 whose columns are the values of
      %  the Hermite base.
      %
      %  .. math::
      %
      %    \begin{eqnarray}
      %    h_1(t) &=& x^2(3-2x)  \\
      %    h_0(t) &=& 1-h_1(t)    \\
      %    h_2(t) &=& x(x(x-2)+1) \\
      %    h_3(t) &=& x^2(x-1)
      %    \end{eqnarray}
      %
      %  basis can be returned in separated vectors
      %
      %  .. code-block:: matlab
      %
      %     [h0,h1,h2,h3] = object.base( t )      % base for the interval [0,1]
      %     [h0,h1,h2,h3] = object.base( t, H )   % base for the interval [0,H]
      %
      [varargout{1:nargout}] = BaseHermiteWrapper('base',varargin{:});
    end
    % --------------------------------------------------------------------
    function varargout = base_D( ~, varargin )1
      %
      %  Evaluate an Hermite base derivative at point(s) t for the interval 
      %  :math:`[0,1]` (or :math:`[0,H]` if second argument
      %  is present).
      %
      %  .. code-block:: matlab
      %
      %      BASE = object.base_D( t )      % base derivative for the interval [0,1]
      %      BASE = object.base_D( t, H )   % base derivative for the interval [0,H]
      %
      %  BASE is a matrix length(t) x 4 whose columns are the values of
      %  the Hermite base.
      %
      %  basis can be returned in separated vectors
      %
      %  .. code-block:: matlab
      %
      %     [h0,h1,h2,h3] = object.base_D( t )      % base derivative for the interval [0,1]
      %     [h0,h1,h2,h3] = object.base_D( t, H )   % base derivative for the interval [0,H]
      %
      [varargout{1:nargout}] = BaseHermiteWrapper('base_D',varargin{:});
    end
    % --------------------------------------------------------------------
    function varargout = base_DD( ~, varargin )
      %
      %  Evaluate an Hermite base second derivative at point(s) t for the interval 
      %  :math:`[0,1]` (or :math:`[0,H]` if second argument
      %  is present).
      %
      %  .. code-block:: matlab
      %
      %      BASE = object.base_DD( t )      % base second derivative for the interval [0,1]
      %      BASE = object.base_DD( t, H )   % base second derivative for the interval [0,H]
      %
      %  BASE is a matrix length(t) x 4 whose columns are the values of
      %  the Hermite base.
      %
      %  basis can be returned in separated vectors
      %
      %  .. code-block:: matlab
      %
      %     [h0,h1,h2,h3] = object.base_DD( t )      % base second derivative for the interval [0,1]
      %     [h0,h1,h2,h3] = object.base_DD( t, H )   % base second derivative for the interval [0,H]
      %
      [varargout{1:nargout}] = BaseHermiteWrapper('base_DD',varargin{:});
    end
    % --------------------------------------------------------------------
    function varargout = base_DDD( ~, varargin )
      %
      %  Evaluate an Hermite base third derivative at point(s) t for the interval 
      %  :math:`[0,1]` (or :math:`[0,H]` if second argument
      %  is present).
      %
      %  .. code-block:: matlab
      %
      %      BASE = object.base_DDD( t )      % base third derivative for the interval [0,1]
      %      BASE = object.base_DDD( t, H )   % base third derivative for the interval [0,H]
      %
      %  BASE is a matrix length(t) x 4 whose columns are the values of
      %  the Hermite base.
      %
      %  basis can be returned in separated vectors
      %
      %  .. code-block:: matlab
      %
      %     [h0,h1,h2,h3] = object.base_DDD( t )      % base third derivative for the interval [0,1]
      %     [h0,h1,h2,h3] = object.base_DDD( t, H )   % base third derivative for the interval [0,H]
      %
      [varargout{1:nargout}] = BaseHermiteWrapper('base_DDD',varargin{:});
    end
    %
    % --------------------------------------------------------------------
    %
    function P = eval( ~, varargin )
      %
      %  Evaluate the cubic polynomial defined on Hermite data:
      %
      %  .. math::
      %
      %   \begin{eqnarray}
      %     \mathbf{p}(t)   &=& h_0(t)\mathbf{p}_0+ h_1(t)\mathbf{p}_1+ h_2(t)\mathbf{t}_0+ h_3(t)\mathbf{t}_1 \\
      %     \mathbf{p}(t,H) &=& h_0(t/H)\mathbf{p}_0+ h_1(t/H)\mathbf{p}_1+ H (h_2(t/H)\mathbf{t}_0+ h_3(t/H)\mathbf{t}_1)
      %   \end{eqnarray}
      %
      %  .. code-block:: matlab
      %
      %      values = object.eval( t, P0, P1, T0, T1 )
      %      values = object.eval( t, P0, P1, T0, T1, H )
      %
      P = BaseHermiteWrapper('eval',varargin{:});
    end
    % --------------------------------------------------------------------
    function dP = eval_D( ~, varargin )
      %
      %  Evaluate the derivative :math:`\mathbf{p}'(t)` of the cubic polynomial defined on Hermite data:
      %
      %  .. code-block:: matlab
      %
      %      values = object.eval_D( t, P0, P1, T0, T1 )
      %      values = object.eval_D( t, P0, P1, T0, T1, H )
      %
      dP = BaseHermiteWrapper('eval_D',varargin{:});
    end
    % --------------------------------------------------------------------
    function ddP = eval_DD( ~, varargin )
      %
      %  Evaluate the second derivative :math:`\mathbf{p}''(t)` of the cubic polynomial defined on Hermite data:
      %
      %  .. code-block:: matlab
      %
      %      values = object.eval_DD( t, P0, P1, T0, T1 )
      %      values = object.eval_DD( t, P0, P1, T0, T1, H )
      %
      ddP = BaseHermiteWrapper('eval_DD',varargin{:});
    end
    % --------------------------------------------------------------------
    function dddP = eval_DDD( ~, varargin )
      %
      %  Evaluate the third derivative :math:`\mathbf{p}'''(t)` of the cubic polynomial defined on Hermite data:
      %
      %  .. code-block:: matlab
      %
      %      values = object.eval_DDD( t, P0, P1, T0, T1 )
      %      values = object.eval_DDD( t, P0, P1, T0, T1, H )
      %
      dddP = BaseHermiteWrapper('eval_DDD',varargin{:});
    end
    %
    % --------------------------------------------------------------------
    %
    function varargout = base5( ~, varargin )
      %
      %  Evaluate an Hermite base (quintic degree) at point(s) t for the interval 
      %  :math:`[0,1]` (or :math:`[0,H]` if second argument
      %  is present).
      %
      %  .. code-block:: matlab
      %
      %      BASE = object.base5( t )      % base for the interval [0,1]
      %      BASE = object.base5( t, H )   % base for the interval [0,H]
      %
      %  BASE is a matrix length(t) x 6 whose columns are the values of
      %  the Hermite base that are 6 polynomials of degree 5 with the properties
      %
      %  .. math::
      %
      %    \begin{equation}
      %      \begin{array}{cccccc}
      %        h_1(0) = 1 & h_1(1) = 0 & h'_1(0) = 1 & h'_1(1) = 0 & h''_1(1) = 1 & h''_1(1) = 0 \\
      %        h_0(0) = 0 & h_0(1) = 1 & h'_0(0) = 0 & h'_0(1) = 0 & h''_0(1) = 0 & h''_0(1) = 0 \\
      %        h_2(0) = 0 & h_2(1) = 0 & h'_2(0) = 1 & h'_2(1) = 0 & h''_2(1) = 0 & h''_2(1) = 0 \\
      %        h_3(0) = 0 & h_3(1) = 0 & h'_3(0) = 0 & h'_3(1) = 1 & h''_3(1) = 0 & h''_3(1) = 0 \\
      %        h_4(0) = 0 & h_4(1) = 0 & h'_4(0) = 0 & h'_4(1) = 0 & h''_4(1) = 1 & h''_4(1) = 0 \\
      %        h_5(0) = 0 & h_5(1) = 0 & h'_5(0) = 0 & h'_5(1) = 0 & h''_5(1) = 0 & h''_5(1) = 1 \\
      %      \end{array}
      %    \end{equation}
      %
      %  basis can be returned in separated vectors
      %
      %  .. code-block:: matlab
      %
      %     [h0,h1,h2,h3,h4,h5] = object.base5( t )      % base for the interval [0,1]
      %     [h0,h1,h2,h3,h4,h5] = object.base5( t, H )   % base for the interval [0,H]
      %
      [varargout{1:nargout}] = BaseHermiteWrapper('base5',varargin{:});
    end
    % --------------------------------------------------------------------
    function varargout = base5_D( ~, varargin )
      %
      %  Evaluate an Hermite base derivatives (quintic degree) at point(s) t for the interval 
      %  :math:`[0,1]` (or :math:`[0,H]` if second argument
      %  is present).
      %
      %  .. code-block:: matlab
      %
      %      BASE = object.base5_D( t )      % base for the interval [0,1]
      %      BASE = object.base5_D( t, H )   % base for the interval [0,H]
      %
      %  BASE is a matrix length(t) x 6 whose columns are the values of
      %  the Hermite base.
      %
      %  basis can be returned in separated vectors
      %
      %  .. code-block:: matlab
      %
      %     [h0,h1,h2,h3,h4,h5] = object.base5_D( t )      % base for the interval [0,1]
      %     [h0,h1,h2,h3,h4,h5] = object.base5_D( t, H )   % base for the interval [0,H]
      %
      [varargout{1:nargout}] = BaseHermiteWrapper('base5_D',varargin{:});
    end
    % --------------------------------------------------------------------
    function varargout = base5_DD( ~, varargin )
      %
      %  Evaluate an Hermite base second derivatives (quintic degree) at point(s) t for the interval 
      %  :math:`[0,1]` (or :math:`[0,H]` if second argument
      %  is present).
      %
      %  .. code-block:: matlab
      %
      %      BASE = object.base5_DD( t )      % base for the interval [0,1]
      %      BASE = object.base5_DD( t, H )   % base for the interval [0,H]
      %
      %  BASE is a matrix length(t) x 6 whose columns are the values of
      %  the Hermite base.
      %
      %  basis can be returned in separated vectors
      %
      %  .. code-block:: matlab
      %
      %     [h0,h1,h2,h3,h4,h5] = object.base5_DD( t )      % base for the interval [0,1]
      %     [h0,h1,h2,h3,h4,h5] = object.base5_DD( t, H )   % base for the interval [0,H]
      %
      [varargout{1:nargout}] = BaseHermiteWrapper('base5_DD',varargin{:});
    end
    % --------------------------------------------------------------------
    function varargout = base5_DDD( ~, varargin )
      %
      %  Evaluate an Hermite base third derivatives (quintic degree) at point(s) t for the interval 
      %  :math:`[0,1]` (or :math:`[0,H]` if second argument
      %  is present).
      %
      %  .. code-block:: matlab
      %
      %      BASE = object.base5_DDD( t )      % base for the interval [0,1]
      %      BASE = object.base5_DDD( t, H )   % base for the interval [0,H]
      %
      %  BASE is a matrix length(t) x 6 whose columns are the values of
      %  the Hermite base.
      %
      %  basis can be returned in separated vectors
      %
      %  .. code-block:: matlab
      %
      %     [h0,h1,h2,h3,h4,h5] = object.base5_DDD( t )      % base for the interval [0,1]
      %     [h0,h1,h2,h3,h4,h5] = object.base5_DDD( t, H )   % base for the interval [0,H]
      %
      [varargout{1:nargout}] = BaseHermiteWrapper('base5_DDD',varargin{:});
    end
    % --------------------------------------------------------------------
    function varargout = base5_DDDD( ~, varargin )
      %
      %  Evaluate an Hermite base 4th derivatives (quintic degree) at point(s) t for the interval 
      %  :math:`[0,1]` (or :math:`[0,H]` if second argument
      %  is present).
      %
      %  .. code-block:: matlab
      %
      %      BASE = object.base5_DDDD( t )      % base for the interval [0,1]
      %      BASE = object.base5_DDDD( t, H )   % base for the interval [0,H]
      %
      %  BASE is a matrix length(t) x 6 whose columns are the values of
      %  the Hermite base.
      %
      %  basis can be returned in separated vectors
      %
      %  .. code-block:: matlab
      %
      %     [h0,h1,h2,h3,h4,h5] = object.base5_DDDD( t )      % base for the interval [0,1]
      %     [h0,h1,h2,h3,h4,h5] = object.base5_DDDD( t, H )   % base for the interval [0,H]
      %
      [varargout{1:nargout}] = BaseHermiteWrapper('base5_DDDD',varargin{:});
    end
    % --------------------------------------------------------------------
    function varargout = base5_DDDDD( ~, varargin )
      %
      %  Evaluate an Hermite base 5th derivatives (quintic degree) at point(s) t for the interval 
      %  :math:`[0,1]` (or :math:`[0,H]` if second argument
      %  is present).
      %
      %  .. code-block:: matlab
      %
      %      BASE = object.base5_DDDDD( t )      % base for the interval [0,1]
      %      BASE = object.base5_DDDDD( t, H )   % base for the interval [0,H]
      %
      %  BASE is a matrix length(t) x 6 whose columns are the values of
      %  the Hermite base.
      %
      %  basis can be returned in separated vectors
      %
      %  .. code-block:: matlab
      %
      %     [h0,h1,h2,h3,h4,h5] = object.base5_DDDDD( t )      % base for the interval [0,1]
      %     [h0,h1,h2,h3,h4,h5] = object.base5_DDDDD( t, H )   % base for the interval [0,H]
      %
      [varargout{1:nargout}] = BaseHermiteWrapper('base5_DDDDD',varargin{:});
    end
    %
    % --------------------------------------------------------------------
    %
    function P = eval5( ~, varargin )
      %
      %  Evaluate the quintic polynomial defined on Hermite data:
      %
      %  .. math::
      %
      %   \begin{eqnarray}
      %     \mathbf{p}(t)   &=& h_0(t)\mathbf{p}_0+ h_1(t)\mathbf{p}_1+ 
      %                         h_2(t)\mathbf{t}_0+ h_3(t)\mathbf{t}_1+
      %                         h_4(t)\mathbf{j}_0+ h_5(t)\mathbf{j}_1  \\
      %     \mathbf{p}(t,H) &=& h_0(t/H)\mathbf{p}_0+ h_1(t/H)\mathbf{p}_1+
      %                         H (h_2(t/H)\mathbf{t}_0+ h_3(t/H)\mathbf{t}_1)+
      %                         H^2 (h_4(t/H)\mathbf{t}_0+ h_5(t/H)\mathbf{j}_1)
      %   \end{eqnarray}
      %
      %  .. code-block:: matlab
      %
      %      values = object.eval5( t, P0, P1, T0, T1, J0, J1 )
      %      values = object.eval5( t, P0, P1, T0, T1, J0, J1, H )
      %
      P = BaseHermiteWrapper('eval5',varargin{:});
    end
    % --------------------------------------------------------------------
    function dP = eval5_D( ~, varargin )
      %
      %  Evaluate the derivative :math:`\mathbf{p}'(t)` of the quintic polynomial defined on Hermite data:
      %
      %  .. code-block:: matlab
      %
      %      values = object.eval5_D( t, P0, P1, T0, T1, J0, J1 )
      %      values = object.eval5_D( t, P0, P1, T0, T1, J0, J1, H )
      %
      dP = BaseHermiteWrapper('eval5_D',varargin{:});
    end
    % --------------------------------------------------------------------
    function ddP = eval5_DD( ~, varargin )
      %
      %  Evaluate the second derivative :math:`\mathbf{p}''(t)` of the quintic polynomial defined on Hermite data:
      %
      %  .. code-block:: matlab
      %
      %      values = object.eval5_DD( t, P0, P1, T0, T1, J0, J1 )
      %      values = object.eval5_DD( t, P0, P1, T0, T1, J0, J1, H )
      %
      ddP = BaseHermiteWrapper('eval5_DD',varargin{:});
    end
    % --------------------------------------------------------------------
    function dddP = eval5_DDD( ~, varargin )
      %
      %  Evaluate the third derivative :math:`\mathbf{p}'''(t)` of the quintic polynomial defined on Hermite data:
      %
      %  .. code-block:: matlab
      %
      %      values = object.eval5_DDD( t, P0, P1, T0, T1, J0, J1 )
      %      values = object.eval5_DDD( t, P0, P1, T0, T1, J0, J1, H )
      %
      dddP = BaseHermiteWrapper('eval5_DDD',varargin{:});
    end
    % --------------------------------------------------------------------
    function ddddP = eval5_DDDD( ~, varargin )
      %
      %  Evaluate the 4th derivative :math:`\mathbf{p}''''(t)` of the quintic polynomial defined on Hermite data:
      %
      %  .. code-block:: matlab
      %
      %      values = object.eval5_DDDD( t, P0, P1, T0, T1, J0, J1 )
      %      values = object.eval5_DDDD( t, P0, P1, T0, T1, J0, J1, H )
      %
      ddddP = BaseHermiteWrapper('eval5_DDDD',varargin{:});
    end
    % --------------------------------------------------------------------
    function dddddP = eval5_DDDDD( ~, varargin )
      %
      %  Evaluate the 5th derivative :math:`\mathbf{p}'''''(t)` of the quintic polynomial defined on Hermite data:
      %
      %  .. code-block:: matlab
      %
      %      values = object.eval5_DDDDD( t, P0, P1, T0, T1, J0, J1 )
      %      values = object.eval5_DDDDD( t, P0, P1, T0, T1, J0, J1, H )
      %
      dddddP = BaseHermiteWrapper('eval5_DDDDD',varargin{:});
    end
    %
    % --------------------------------------------------------------------
    %
    function [P0,P1,P2,P3] = hermite_to_bezier( ~, p0, p1, t0, t1 )
      %
      % Given the cubic polynomial defined on Hermite data:
      %
      %  .. math::
      %
      %     \mathbf{p}(t) = h_0(t)\mathbf{p}_0+ h_1(t)\mathbf{p}_1+ h_2(t)\mathbf{t}_0+ h_3(t)\mathbf{t}_1
      %
      %  return the Bezier polynomial of the same cubic
      %
      [P0,P1,P2,P3] = BaseHermiteWrapper('hermite_to_bezier',p0, p1, t0, t1);
    end
    % --------------------------------------------------------------------
    function [P0,P1,T0,T1] = bezier_to_hermite( ~, p0, p1, p2, p3 )
      %
      % Given the cubic polynomial defined with Bezier polygon
      % return the Hermite data for the same polynomial
      %
      [P0,P1,T0,T1] = BaseHermiteWrapper('bezier_to_hermite',p0, p1, p2, p3);
    end
    %
    % --------------------------------------------------------------------
    %
    function [D1,sqrtD1] = L2_first_derivative( ~ )
      sqrtD1 = BaseHermiteWrapper('L2_first_derivative');
      D1     = sqrtD1*sqrtD1.';
    end
    % --------------------------------------------------------------------
    function [D2,sqrtD2] = L2_second_derivative( ~ )
      sqrtD2 = BaseHermiteWrapper('L2_second_derivative');
      D2     = sqrtD2*sqrtD2.';
    end
    % --------------------------------------------------------------------
    function [D3,sqrtD3] = L2_third_derivative( ~ )
      sqrtD3 = BaseHermiteWrapper('L2_third_derivative');
      D3     = sqrtD3*sqrtD3.';
    end
    % --------------------------------------------------------------------
    function L = approximate_length( ~, varargin )
      %
      %  Approximate the length of the cubic polynomial :math:`\mathbf{p}(t)` defined on Hermite data:
      %
      %  .. code-block:: matlab
      %
      %      values = object.approximate_length( t, P0, P1, T0, T1 )
      %      values = object.approximate_length( t, P0, P1, T0, T1, H )
      %
      %  The length is approximated usin 100 linear segment.
      %
      L = BaseHermiteWrapper( 'approximate_length', varargin{:} );
    end
    % --------------------------------------------------------------------
    function [P0,P1,T0,T1] = cut( ~, varargin )
      %
      %  Cut the cubic polynomial :math:`\mathbf{p}(t)` defined on Hermite data
      %  on the interval [a,b] and return the new Hermite data
      %
      %  .. code-block:: matlab
      %
      %      [new_P0,new_P1,new_T0,new_T1] = object.cut( a, b, P0, P1, T0, T1 )
      %      [new_P0,new_P1,new_T0,new_T1] = object.cut( a, b, P0, P1, T0, T1, H )
      %
      %  The parametrization of the new Hermite data is on [0,1]
      %
      [P0,P1,T0,T1] = BaseHermiteWrapper( 'cut', varargin{:} );
    end
  end
end
