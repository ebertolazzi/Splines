classdef BaseHermite < handle

  methods
    %>
    %> Build a matlab object storing Hermite base evaluator
    %>
    function self = BaseHermite()
    end

    function delete(self)
    end
    % --------------------------------------------------------------------
    %>
    %> Evaluate an Hermite base (cubic degree) at point(s) t for the interval
    %> \f$ [0,1] \f$ (or \f$ [0,H] \f$ if second argument
    %> is present).
    %>
    %> ```{matlab}
    %>   BASE = object.base( t )      % base for the interval [0,1]
    %>   BASE = object.base( t, H )   % base for the interval [0,H]
    %> ```
    %>
    %> BASE is a matrix length(t) x 4 whose columns are the values of
    %> the Hermite base.
    %>
    %>
    %> \f{eqnarray*}{
    %>     h_1(t) &=& x^2(3-2x)  \\
    %>     h_0(t) &=& 1-h_1(t)    \\
    %>     h_2(t) &=& x(x(x-2)+1) \\
    %>     h_3(t) &=& x^2(x-1)
    %> \f}
    %>
    %> basis can be returned in separated vectors
    %>
    %> ```{matlab}
    %>   [h0,h1,h2,h3] = object.base( t )      % base for the interval [0,1]
    %>   [h0,h1,h2,h3] = object.base( t, H )   % base for the interval [0,H]
    %> ```
    %>
    function varargout = base( ~, varargin )
      [varargout{1:nargout}] = BaseHermiteMexWrapper('base',varargin{:});
    end
    % --------------------------------------------------------------------
    %>
    %> Evaluate an Hermite base derivative at point(s) t for the interval
    %> \f$ [0,1] \f$ (or \f$ [0,H] \f$ if second argument
    %> is present).
    %>
    %> ```{matlab}
    %>   BASE = object.base_D( t )      % base derivative for the interval [0,1]
    %>   BASE = object.base_D( t, H )   % base derivative for the interval [0,H]
    %> ```
    %>
    %> BASE is a matrix length(t) x 4 whose columns are the values of
    %> the Hermite base.
    %>
    %> basis can be returned in separated vectors
    %>
    %> ```{matlab}
    %>   [h0,h1,h2,h3] = object.base_D( t )      % base derivative for the interval [0,1]
    %>   [h0,h1,h2,h3] = object.base_D( t, H )   % base derivative for the interval [0,H]
    %> ```
    %>
    function varargout = base_D( ~, varargin )
      [varargout{1:nargout}] = BaseHermiteMexWrapper('base_D',varargin{:});
    end
    % --------------------------------------------------------------------
    %>
    %> Evaluate an Hermite base second derivative at point(s) t for the interval
    %> \f$ [0,1] \f$ (or \f$ [0,H] \f$ if second argument is present).
    %>
    %> ```{matlab}
    %>   BASE = object.base_DD( t )      % base second derivative for the interval [0,1]
    %>   BASE = object.base_DD( t, H )   % base second derivative for the interval [0,H]
    %> ```
    %>
    %> BASE is a matrix length(t) x 4 whose columns are the values of
    %> the Hermite base.
    %>
    %> basis can be returned in separated vectors
    %>
    %> ```{matlab}
    %>   [h0,h1,h2,h3] = object.base_DD( t )      % base second derivative for the interval [0,1]
    %>   [h0,h1,h2,h3] = object.base_DD( t, H )   % base second derivative for the interval [0,H]
    %> ```
    %>
    function varargout = base_DD( ~, varargin )
      [varargout{1:nargout}] = BaseHermiteMexWrapper('base_DD',varargin{:});
    end
    % --------------------------------------------------------------------
    %>
    %> Evaluate an Hermite base third derivative at point(s) t for the interval
    %> \f$ [0,1] \f$ (or \f$ [0,H] \f$ if second argument is present).
    %>
    %> ```{matlab}
    %>   BASE = object.base_DDD( t )      % base third derivative for the interval [0,1]
    %>   BASE = object.base_DDD( t, H )   % base third derivative for the interval [0,H]
    %> ```
    %>
    %> BASE is a matrix length(t) x 4 whose columns are the values of the Hermite base.
    %>
    %> basis can be returned in separated vectors
    %>
    %> ```{matlab}
    %>   [h0,h1,h2,h3] = object.base_DDD( t )      % base third derivative for the interval [0,1]
    %>   [h0,h1,h2,h3] = object.base_DDD( t, H )   % base third derivative for the interval [0,H]
    %> ```
    %>
    function varargout = base_DDD( ~, varargin )
      [varargout{1:nargout}] = BaseHermiteMexWrapper('base_DDD',varargin{:});
    end
    % --------------------------------------------------------------------
    %>
    %> Evaluate the cubic polynomial defined on Hermite data:
    %>
    %>
    %> \f{eqnarray*}{
    %>   \mathbf{p}(t)   &=& h_0(t)\mathbf{p}_0+ h_1(t)\mathbf{p}_1+ h_2(t)\mathbf{t}_0+ h_3(t)\mathbf{t}_1 \\
    %>   \mathbf{p}(t,H) &=& h_0(t/H)\mathbf{p}_0+ h_1(t/H)\mathbf{p}_1+ H (h_2(t/H)\mathbf{t}_0+ h_3(t/H)\mathbf{t}_1)
    %> \f}
    %>
    %> ```{matlab}
    %>   values = object.eval( t, P0, P1, T0, T1 )
    %>   values = object.eval( t, P0, P1, T0, T1, H )
    %> ```
    %>
    function P = eval( ~, varargin )
      P = BaseHermiteMexWrapper('eval',varargin{:});
    end
    % --------------------------------------------------------------------
    %>
    %> Evaluate the derivative \f$ \mathbf{p}'(t) \f$ of the cubic
    %> polynomial defined on Hermite data:
    %>
    %>
    %> ```{matlab}
    %>   values = object.eval_D( t, P0, P1, T0, T1 )
    %>   values = object.eval_D( t, P0, P1, T0, T1, H )
    %> ```
    %>
    function dP = eval_D( ~, varargin )
      dP = BaseHermiteMexWrapper('eval_D',varargin{:});
    end
    % --------------------------------------------------------------------
    %>
    %> Evaluate the second derivative \f$ \mathbf{p}''(t) \f$ of
    %> the cubic polynomial defined on Hermite data:
    %>
    %> ```{matlab}
    %>   values = object.eval_DD( t, P0, P1, T0, T1 )
    %>   values = object.eval_DD( t, P0, P1, T0, T1, H )
    %> ```
    %>
    function ddP = eval_DD( ~, varargin )
      ddP = BaseHermiteMexWrapper('eval_DD',varargin{:});
    end
    % --------------------------------------------------------------------
    %>
    %> Evaluate the third derivative \f$ \mathbf{p}'''(t) \f$ of the
    %> cubic polynomial defined on Hermite data:
    %>
    %> ```{matlab}
    %>   values = object.eval_DDD( t, P0, P1, T0, T1 )
    %>   values = object.eval_DDD( t, P0, P1, T0, T1, H )
    %> ```
    %>
    function dddP = eval_DDD( ~, varargin )
      dddP = BaseHermiteMexWrapper('eval_DDD',varargin{:});
    end
    % --------------------------------------------------------------------
    %>
    %> Evaluate an Hermite base (quintic degree) at point(s) t for the interval
    %> \f$ [0,1] \f$ (or \f$ [0,H] \f$ if second argument is present).
    %>
    %> ```{matlab}
    %>   BASE = object.base5( t )      % base for the interval [0,1]
    %>   BASE = object.base5( t, H )   % base for the interval [0,H]
    %> ```
    %>
    %> BASE is a matrix length(t) x 6 whose columns are the values of
    %> the Hermite base that are 6 polynomials of degree 5 with the properties
    %>
    %>
    %> \f[
    %>       \begin{array}{cccccc}
    %>         h_1(0) = 1 & h_1(1) = 0 & h'_1(0) = 1 & h'_1(1) = 0 & h''_1(1) = 1 & h''_1(1) = 0 \\
    %>         h_0(0) = 0 & h_0(1) = 1 & h'_0(0) = 0 & h'_0(1) = 0 & h''_0(1) = 0 & h''_0(1) = 0 \\
    %>         h_2(0) = 0 & h_2(1) = 0 & h'_2(0) = 1 & h'_2(1) = 0 & h''_2(1) = 0 & h''_2(1) = 0 \\
    %>         h_3(0) = 0 & h_3(1) = 0 & h'_3(0) = 0 & h'_3(1) = 1 & h''_3(1) = 0 & h''_3(1) = 0 \\
    %>         h_4(0) = 0 & h_4(1) = 0 & h'_4(0) = 0 & h'_4(1) = 0 & h''_4(1) = 1 & h''_4(1) = 0 \\
    %>         h_5(0) = 0 & h_5(1) = 0 & h'_5(0) = 0 & h'_5(1) = 0 & h''_5(1) = 0 & h''_5(1) = 1 \\
    %>       \end{array}
    %> \f]
    %>
    %> basis can be returned in separated vectors
    %>
    %> ```{matlab}
    %>   [h0,h1,h2,h3,h4,h5] = object.base5( t )      % base for the interval [0,1]
    %>   [h0,h1,h2,h3,h4,h5] = object.base5( t, H )   % base for the interval [0,H]
    %> ```
    %>
    function varargout = base5( ~, varargin )
      [varargout{1:nargout}] = BaseHermiteMexWrapper('base5',varargin{:});
    end
    % --------------------------------------------------------------------
    %>
    %> Evaluate an Hermite base derivatives (quintic degree) at point(s) t for the interval
    %> \f$ [0,1] \f$ (or \f$ [0,H] \f$ if second argument is present).
    %>
    %> ```{matlab}
    %>   BASE = object.base5_D( t )      % base for the interval [0,1]
    %>   BASE = object.base5_D( t, H )   % base for the interval [0,H]
    %> ```
    %>
    %> BASE is a matrix length(t) x 6 whose columns are the values of
    %> the Hermite base.
    %>
    %> basis can be returned in separated vectors
    %>
    %> ```{matlab}
    %>   [h0,h1,h2,h3,h4,h5] = object.base5_D( t )      % base for the interval [0,1]
    %>   [h0,h1,h2,h3,h4,h5] = object.base5_D( t, H )   % base for the interval [0,H]
    %> ```
    %>
    function varargout = base5_D( ~, varargin )
      [varargout{1:nargout}] = BaseHermiteMexWrapper('base5_D',varargin{:});
    end
    % --------------------------------------------------------------------
    %>
    %> Evaluate an Hermite base second derivatives (quintic degree) at point(s) t for the interval
    %> \f$ [0,1] \f$ (or \f$ [0,H] \f$ if second argument is present).
    %>
    %>
    %> ```{matlab}
    %>   BASE = object.base5_DD( t )      % base for the interval [0,1]
    %>   BASE = object.base5_DD( t, H )   % base for the interval [0,H]
    %> ```
    %>
    %> BASE is a matrix length(t) x 6 whose columns are the values of
    %> the Hermite base.
    %>
    %> basis can be returned in separated vectors
    %>
    %> ```{matlab}
    %>   [h0,h1,h2,h3,h4,h5] = object.base5_DD( t )      % base for the interval [0,1]
    %>   [h0,h1,h2,h3,h4,h5] = object.base5_DD( t, H )   % base for the interval [0,H]
    %> ```
    %>
    function varargout = base5_DD( ~, varargin )
      [varargout{1:nargout}] = BaseHermiteMexWrapper('base5_DD',varargin{:});
    end
    % --------------------------------------------------------------------
    %>
    %> Evaluate an Hermite base third derivatives (quintic degree) at point(s)
    %> t for the interval \f$ [0,1] \f$ (or \f$ [0,H] \f$ if second argument
    %> is present).
    %>
    %> ```{matlab}
    %>   BASE = object.base5_DDD( t )      % base for the interval [0,1]
    %>   BASE = object.base5_DDD( t, H )   % base for the interval [0,H]
    %> ```
    %>
    %> BASE is a matrix length(t) x 6 whose columns are the values of
    %> the Hermite base.
    %>
    %> basis can be returned in separated vectors
    %>
    %> ```{matlab}
    %>  [h0,h1,h2,h3,h4,h5] = object.base5_DDD( t )      % base for the interval [0,1]
    %>  [h0,h1,h2,h3,h4,h5] = object.base5_DDD( t, H )   % base for the interval [0,H]
    %> ```
    %>
    function varargout = base5_DDD( ~, varargin )
      [varargout{1:nargout}] = BaseHermiteMexWrapper('base5_DDD',varargin{:});
    end
    % --------------------------------------------------------------------
    %>
    %> Evaluate an Hermite base 4th derivatives (quintic degree) at point(s)
    %> t for the interval \f$ [0,1] \f$ (or \f$ [0,H] \f$ if second argument
    %> is present).
    %>
    %> ```{matlab}
    %>   BASE = object.base5_DDDD( t )      % base for the interval [0,1]
    %>   BASE = object.base5_DDDD( t, H )   % base for the interval [0,H]
    %> ```
    %>
    %> BASE is a matrix length(t) x 6 whose columns are the values of
    %> the Hermite base.
    %>
    %> basis can be returned in separated vectors
    %>
    %> ```{matlab}
    %>   [h0,h1,h2,h3,h4,h5] = object.base5_DDDD( t )      % base for the interval [0,1]
    %>   [h0,h1,h2,h3,h4,h5] = object.base5_DDDD( t, H )   % base for the interval [0,H]
    %> ```
    %>
    function varargout = base5_DDDD( ~, varargin )
      [varargout{1:nargout}] = BaseHermiteMexWrapper('base5_DDDD',varargin{:});
    end
    % --------------------------------------------------------------------
    %>
    %> Evaluate an Hermite base 5th derivatives (quintic degree) at point(s)
    %> t for the interval \f$ [0,1] \f$ (or \f$ [0,H] \f$ if second argument
    %> is present).
    %>
    %> ```{matlab}
    %>   BASE = object.base5_DDDDD( t )      % base for the interval [0,1]
    %>   BASE = object.base5_DDDDD( t, H )   % base for the interval [0,H]
    %> ```
    %>
    %> BASE is a matrix length(t) x 6 whose columns are the values of
    %> the Hermite base.
    %>
    %> basis can be returned in separated vectors
    %>
    %> ```{matlab}
    %>   [h0,h1,h2,h3,h4,h5] = object.base5_DDDDD( t )      % base for the interval [0,1]
    %>   [h0,h1,h2,h3,h4,h5] = object.base5_DDDDD( t, H )   % base for the interval [0,H]
    %> ```
    %>
    function varargout = base5_DDDDD( ~, varargin )
      [varargout{1:nargout}] = BaseHermiteMexWrapper('base5_DDDDD',varargin{:});
    end
    % --------------------------------------------------------------------
    %>
    %> Evaluate the quintic polynomial defined on Hermite data:
    %>
    %> \f{eqnarray*}{
    %>       \mathbf{p}(t)   &=& h_0(t)\mathbf{p}_0+ h_1(t)\mathbf{p}_1+
    %>                           h_2(t)\mathbf{t}_0+ h_3(t)\mathbf{t}_1+
    %>                           h_4(t)\mathbf{j}_0+ h_5(t)\mathbf{j}_1  \\
    %>       \mathbf{p}(t,H) &=& h_0(t/H)\mathbf{p}_0+ h_1(t/H)\mathbf{p}_1+
    %>                           H (h_2(t/H)\mathbf{t}_0+ h_3(t/H)\mathbf{t}_1)+
    %>                           H^2 (h_4(t/H)\mathbf{t}_0+ h_5(t/H)\mathbf{j}_1)
    %> \f}
    %>
    %> ```{matlab}
    %>   values = object.eval5( t, P0, P1, T0, T1, J0, J1 )
    %>   values = object.eval5( t, P0, P1, T0, T1, J0, J1, H )
    %> ```
    %>
    function P = eval5( ~, varargin )
      P = BaseHermiteMexWrapper('eval5',varargin{:});
    end
    % --------------------------------------------------------------------
    %>
    %> Evaluate the derivative \f$ \mathbf{p}'(t) \f$ of the
    %> quintic polynomial defined on Hermite data:
    %>
    %> ```{matlab}
    %>   values = object.eval5_D( t, P0, P1, T0, T1, J0, J1 )
    %>   values = object.eval5_D( t, P0, P1, T0, T1, J0, J1, H )
    %> ```
    %>
    function dP = eval5_D( ~, varargin )
      dP = BaseHermiteMexWrapper('eval5_D',varargin{:});
    end
    % --------------------------------------------------------------------
    %>
    %> Evaluate the second derivative \f$ \mathbf{p}''(t) \f$ of
    %> the quintic polynomial defined on Hermite data:
    %>
    %> ```{matlab}
    %>   values = object.eval5_DD( t, P0, P1, T0, T1, J0, J1 )
    %>   values = object.eval5_DD( t, P0, P1, T0, T1, J0, J1, H )
    %> ```
    %>
    function ddP = eval5_DD( ~, varargin )
      ddP = BaseHermiteMexWrapper('eval5_DD',varargin{:});
    end
    % --------------------------------------------------------------------
    %>
    %> Evaluate the third derivative \f$  \mathbf{p}'''(t) \f$ of
    %> the quintic polynomial defined on Hermite data:
    %>
    %> ```{matlab}
    %>   values = object.eval5_DDD( t, P0, P1, T0, T1, J0, J1 )
    %>   values = object.eval5_DDD( t, P0, P1, T0, T1, J0, J1, H )
    %> ```
    %>
    function dddP = eval5_DDD( ~, varargin )
      dddP = BaseHermiteMexWrapper('eval5_DDD',varargin{:});
    end
    % --------------------------------------------------------------------
    %>
    %> Evaluate the 4th derivative \f$ \mathbf{p}''''(t) \f$ of the
    %> quintic polynomial defined on Hermite data:
    %>
    %> ```{matlab}
    %>   values = object.eval5_DDDD( t, P0, P1, T0, T1, J0, J1 )
    %>   values = object.eval5_DDDD( t, P0, P1, T0, T1, J0, J1, H )
    %> ```
    %>
    function ddddP = eval5_DDDD( ~, varargin )
      ddddP = BaseHermiteMexWrapper('eval5_DDDD',varargin{:});
    end
    % --------------------------------------------------------------------
    %>
    %> Evaluate the 5th derivative \f$ \mathbf{p}'''''(t) \f$ of the
    %> quintic polynomial defined on Hermite data:
    %>
    %> ```{matlab}
    %>   values = object.eval5_DDDDD( t, P0, P1, T0, T1, J0, J1 )
    %>   values = object.eval5_DDDDD( t, P0, P1, T0, T1, J0, J1, H )
    %> ```
    %>
    function dddddP = eval5_DDDDD( ~, varargin )
      dddddP = BaseHermiteMexWrapper('eval5_DDDDD',varargin{:});
    end
    % --------------------------------------------------------------------
    %>
    %> Given the cubic polynomial defined on Hermite data:
    %>
    %> \f{eqnarray*}{
    %>     \mathbf{p}(t) =
    %>     h_0(t)\mathbf{p}_0 + h_1(t)\mathbf{p}_1 +
    %>     h_2(t)\mathbf{t}_0 + h_3(t)\mathbf{t}_1
    %> \f}
    %>
    %> return the Bezier polynomial of the same cubic
    %>
    function [P0,P1,P2,P3] = hermite_to_bezier( ~, p0, p1, t0, t1 )
      [P0,P1,P2,P3] = BaseHermiteMexWrapper('hermite_to_bezier',p0, p1, t0, t1);
    end
    % --------------------------------------------------------------------
    %>
    %> Given the cubic polynomial defined with Bezier polygon
    %> return the Hermite data for the same polynomial
    %>
    function [P0,P1,T0,T1] = bezier_to_hermite( ~, p0, p1, p2, p3 )
      [P0,P1,T0,T1] = BaseHermiteMexWrapper('bezier_to_hermite',p0, p1, p2, p3);
    end
    %
    % --------------------------------------------------------------------
    %
    function [D1,sqrtD1] = L2_first_derivative( ~ )
      sqrtD1 = BaseHermiteMexWrapper('L2_first_derivative');
      D1     = sqrtD1*sqrtD1.';
    end
    % --------------------------------------------------------------------
    function [D2,sqrtD2] = L2_second_derivative( ~ )
      sqrtD2 = BaseHermiteMexWrapper('L2_second_derivative');
      D2     = sqrtD2*sqrtD2.';
    end
    % --------------------------------------------------------------------
    function [D3,sqrtD3] = L2_third_derivative( ~ )
      sqrtD3 = BaseHermiteMexWrapper('L2_third_derivative');
      D3     = sqrtD3*sqrtD3.';
    end
    % --------------------------------------------------------------------
    %>
    %> Approximate the length of the cubic polynomial \f$ \mathbf{p}(t) \f$
    %> defined on Hermite data:
    %>
    %> ```{matlab}
    %>   values = object.approximate_length( t, P0, P1, T0, T1 )
    %>   values = object.approximate_length( t, P0, P1, T0, T1, H )
    %> ```
    %>
    %> The length is approximated usin 100 linear segment.
    %>
    function L = approximate_length( ~, varargin )
      L = BaseHermiteMexWrapper( 'approximate_length', varargin{:} );
    end
    % --------------------------------------------------------------------
    %>
    %> Cut the cubic polynomial \f$ \mathbf{p}(t) \f$ defined on Hermite data
    %> on the interval [a,b] and return the new Hermite data
    %>
    %> ```{matlab}
    %>   [new_P0,new_P1,new_T0,new_T1] = object.cut( a, b, P0, P1, T0, T1 )
    %>   [new_P0,new_P1,new_T0,new_T1] = object.cut( a, b, P0, P1, T0, T1, H )
    %> ```
    %>
    %> The parametrization of the new Hermite data is on [0,1]
    %>
    function [P0,P1,T0,T1] = cut( ~, varargin )
      [P0,P1,T0,T1] = BaseHermiteMexWrapper( 'cut', varargin{:} );
    end
  end
end
