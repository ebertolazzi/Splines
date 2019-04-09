classdef BaseHermite < handle

  methods

    function self = BaseHermite()
    end

    function delete(self)
    end
    %
    % --------------------------------------------------------------------
    %
    function varargout = base( ~, varargin )
      [varargout{1:nargout}] = BaseHermiteWrapper('base',varargin{:});
    end
    % --------------------------------------------------------------------
    function varargout = base_D( ~, varargin )
      [varargout{1:nargout}] = BaseHermiteWrapper('base_D',varargin{:});
    end
    % --------------------------------------------------------------------
    function varargout = base_DD( ~, varargin )
      [varargout{1:nargout}] = BaseHermiteWrapper('base_DD',varargin{:});
    end
    % --------------------------------------------------------------------
    function varargout = base_DDD( ~, varargin )
      [varargout{1:nargout}] = BaseHermiteWrapper('base_DDD',varargin{:});
    end
    %
    % --------------------------------------------------------------------
    %
    function P = eval( ~, varargin )
      P = BaseHermiteWrapper('eval',varargin{:});
    end
    % --------------------------------------------------------------------
    function dP = eval_D( ~, varargin )
      dP = BaseHermiteWrapper('eval_D',varargin{:});
    end
    % --------------------------------------------------------------------
    function ddP = eval_DD( ~, varargin )
      ddP = BaseHermiteWrapper('eval_DD',varargin{:});
    end
    % --------------------------------------------------------------------
    function dddP = eval_DDD( ~, varargin )
      dddP = BaseHermiteWrapper('eval_DDD',varargin{:});
    end
    %
    % --------------------------------------------------------------------
    %
    function varargout = base5( ~, varargin )
      [varargout{1:nargout}] = BaseHermiteWrapper('base5',varargin{:});
    end
    % --------------------------------------------------------------------
    function varargout = base5_D( ~, varargin )
      [varargout{1:nargout}] = BaseHermiteWrapper('base5_D',varargin{:});
    end
    % --------------------------------------------------------------------
    function varargout = base5_DD( ~, varargin )
      [varargout{1:nargout}] = BaseHermiteWrapper('base5_DD',varargin{:});
    end
    % --------------------------------------------------------------------
    function varargout = base5_DDD( ~, varargin )
      [varargout{1:nargout}] = BaseHermiteWrapper('base5_DDD',varargin{:});
    end
    % --------------------------------------------------------------------
    function varargout = base5_DDDD( ~, varargin )
      [varargout{1:nargout}] = BaseHermiteWrapper('base5_DDDD',varargin{:});
    end
    % --------------------------------------------------------------------
    function varargout = base5_DDDDD( ~, varargin )
      [varargout{1:nargout}] = BaseHermiteWrapper('base5_DDDDD',varargin{:});
    end
    %
    % --------------------------------------------------------------------
    %
    function P = eval5( ~, varargin )
      P = BaseHermiteWrapper('eval5',varargin{:});
    end
    % --------------------------------------------------------------------
    function dP = eval5_D( ~, varargin )
      dP = BaseHermiteWrapper('eval5_D',varargin{:});
    end
    % --------------------------------------------------------------------
    function ddP = eval5_DD( ~, varargin )
      ddP = BaseHermiteWrapper('eval5_DD',varargin{:});
    end
    % --------------------------------------------------------------------
    function dddP = eval5_DDD( ~, varargin )
      dddP = BaseHermiteWrapper('eval5_DDD',varargin{:});
    end
    % --------------------------------------------------------------------
    function ddddP = eval5_DDDD( ~, varargin )
      ddddP = BaseHermiteWrapper('eval5_DDDD',varargin{:});
    end
    % --------------------------------------------------------------------
    function dddddP = eval5_DDDDD( ~, varargin )
      dddddP = BaseHermiteWrapper('eval5_DDDDD',varargin{:});
    end
    %
    % --------------------------------------------------------------------
    %
    function [P0,P1,P2,P3] = hermite_to_bezier( ~, p0, p1, t0, t1 )
      [P0,P1,P2,P3] = BaseHermiteWrapper('hermite_to_bezier',p0, p1, t0, t1);
    end
    % --------------------------------------------------------------------
    function [P0,P1,T0,T1] = bezier_to_hermite( ~, p0, p1, p2, p3 )
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
      L = BaseHermiteWrapper( 'approximate_length', varargin{:} );
    end
    % --------------------------------------------------------------------
    function [P0,P1,T0,T1] = cut( ~, varargin )
      [P0,P1,T0,T1] = BaseHermiteWrapper( 'cut', varargin{:} );
    end
  end
end
