function [m_slopes_vec, var_m] = var_est_low(varargin)

if(nargin == 2)
    H_app = varargin{1};
    M = varargin{2};
    dof = 1;
    var_m = zeros(length(H_app),1);
    m_slopes_vec = zeros(length(H_app),1);
    
    for i = dof+1 : length(H_app)
        
        H_matrix = H_app(1:i);
        M_vec = M(1:i);
        m_slopes = H_matrix\M_vec;
        m_slopes_vec(i,:) = m_slopes';
        var_m(i) = 1/(2*i-dof)*sum((H_matrix*m_slopes-M_vec).^2);
    end
    
    
elseif(nargin == 3)
    H_app = varargin{1};
    M_ax = varargin{2};
    M_rad = varargin{3};
    
    dof = 2;
    var_m = zeros(length(H_app),1);
    m_slopes_vec = zeros(length(H_app),2);
    
    for i = dof+1 : length(H_app)
        
        H_matrix = [H_app(1:i) zeros(size(H_app(1:i)));...
            zeros(size(H_app(1:i))) H_app(1:i)];
        M_vec = [M_ax(1:i); M_rad(1:i)];
        m_slopes = H_matrix\M_vec;
        m_slopes_vec(i,:) = m_slopes';
        var_m(i) = 1/(2*i-dof)*sum((H_matrix*m_slopes-M_vec).^2);
    end
    
end

end
