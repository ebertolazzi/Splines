classdef BaseHermite < handle

  methods

    function self = BaseHermite()
    end

    function delete(self)
    end

    function varargout = base( ~, tt, varargin )
      [varargout{1:nargout}] = BaseHermiteWrapper('base',tt,varargin{:});
    end

    function varargout = base_D( ~, tt, varargin )
      [varargout{1:nargout}] = BaseHermiteWrapper('base_D',tt,varargin{:});
    end

    function varargout = base_DD( ~, tt, varargin )
      [varargout{1:nargout}] = BaseHermiteWrapper('base_DD',tt,varargin{:});
    end

    function varargout = base_DDD( ~, tt, varargin )
      [varargout{1:nargout}] = BaseHermiteWrapper('base_DDD',tt,varargin{:});
    end

    function [P0,P1,P2,P3] = hermite_to_bezier( ~, p0, p1, t0, t1 )
      [P0,P1,P2,P3] = BaseHermiteWrapper('hermite_to_bezier',p0, p1, t0, t1);
    end

    function [P0,P1,T0,T1] = bezier_to_hermite( ~, p0, p1, p2, p3 )
      [P0,P1,T0,T1] = BaseHermiteWrapper('bezier_to_hermite',p0, p1, p2, p3);
    end

    function [D1,sqrtD1] = L2_first_derivative( ~ )
      [D1,sqrtD1]= BaseHermiteWrapper('L2_first_derivative');
    end

    function [D2,sqrtD2] = L2_second_derivative( ~ )
      [D2,sqrtD2]= BaseHermiteWrapper('L2_second_derivative');
    end

    function [D3,sqrtD3] = L2_third_derivative( ~ )
      [D2,sqrtD2]= BaseHermiteWrapper('L2_third_derivative');
    end

    function varargout = base5( ~, tt, varargin )
      [varargout{1:nargout}] = BaseHermiteWrapper('base5',tt,varargin{:});
    end

    function varargout = base5_D( ~, tt, varargin )
      [varargout{1:nargout}] = BaseHermiteWrapper('base5_D',tt,varargin{:});
    end

    function varargout = base5_DD( ~, tt, varargin )
      [varargout{1:nargout}] = BaseHermiteWrapper('base5_DD',tt,varargin{:});
    end

    function varargout = base5_DDD( ~, tt, varargin )
      [varargout{1:nargout}] = BaseHermiteWrapper('base5_DDD',tt,varargin{:});
    end

    function varargout = base5_DDDD( ~, tt, varargin )
      [varargout{1:nargout}] = BaseHermiteWrapper('base5_DDDD',tt,varargin{:});
    end

    function varargout = base5_DDDDD( ~, tt, varargin )
      [varargout{1:nargout}] = BaseHermiteWrapper('base5_DDDDD',tt,varargin{:});
    end

  end
end
