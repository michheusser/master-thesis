%function [H_app_simplified, M_simplified ] = SimplifyMH(varargin)

function varargout = SimplifyMH(varargin)

H_app = varargin{1};
M = varargin{2};

if(strcmp(varargin{4},'Zero'))
    m_scale = varargin{3};
elseif(strcmp(varargin{4},'Saturates'))
    M_s = varargin{3};
end

%% Division in 3 Parts


diff_H_app = diff(H_app);

n = 1;
for i = 1 : length(diff_H_app)-1
    if(diff_H_app(i)*diff_H_app(i+1) <= 0)
        I_lim(n) = i+1;
        n = n + 1;
    end
    
end

H_app_start = H_app(1:I_lim(1));
H_app_pos = flipud(H_app(I_lim(1):I_lim(2)));
H_app_neg = H_app(I_lim(2):end);

M_start = M(1:I_lim(1));
M_pos = flipud(M(I_lim(1):I_lim(2)));
M_neg = M(I_lim(2):end);

%% Averaging

H_app_av = H_app_pos;
M_av = 0.5*(M_pos + M_neg);

%% y-Shifting

M_s_pos = max(M_av);
M_s_neg = min(M_av);

M_shift_y = M_av - 0.5*(M_s_pos + M_s_neg);

%% x-Shifting
circshiftM = circshift(M_shift_y,[-1,1]);
M_shift_y_1 = circshiftM(1:end-1);
M_shift_y_0 = M_shift_y(1:end-1);
I = find(M_shift_y_0.*M_shift_y_1<0);


i_cross = I(round(length(I)/2));

for i = 1 : length(M_shift_y)-1
    if(i == i_cross)
        
        m = (M_shift_y(i+1)-M_shift_y(i))/(H_app_av(i+1)-H_app_av(i));
        Dx = -M_shift_y(i)/m + H_app_av(i);
        
    end
    
end

H_app_shift_x = H_app_av - Dx;


%% Scaling to Correct M_s
if(nargin == 3)
    M_s = varargin{3};
    M_scaled = M_shift_y*(M_s/max(M_shift_y));
else
    M_scaled = M_shift_y;
end

    
if(strcmp(varargin{4},'Saturates'))
       m_scale = M_s/max(M_shift_y);
end

M_scaled = M_shift_y*m_scale; 

%% Ideal Curve

H_app_simplified = H_app_shift_x;
M_simplified = M_scaled;

plot(H_app,M,'-*',H_app_simplified,M_simplified,'-*')
grid on
axis tight
legend('Raw','Simplified')

if(strcmp(varargin{4},'Zero'))
    varargout{1} = H_app_simplified;
    varargout{2} = M_simplified;
elseif(strcmp(varargin{4},'Saturates'))
    varargout{1} =  H_app_simplified;
    varargout{2} =  M_simplified;
    varargout{3} =  m_scale;
end

end

