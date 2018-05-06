%% Matched Subspace Detection
% S. Nazari
clear, 
% clc

% Params
SIGNAL  = false;    % Signal in subspace? H1: TRUE ; H0: FALSE
N       = 3;        % Subspace dim
M       = 1000;     % MC runs


f0      = 20;       % Inteference frequency
fss     = 1;
fsc     = 1.4;
phi0    = pi/4;     % Inteference Phase
A       = 2;        % Inteference Amplitude

u       = 2;
a       = 1;
b       = 3;

THRSH   = 6;        % Xi^2 Detection threshold

H1_kn   = 0;
H0_kn   = 0;
H1_uk   = 0;
H0_uk   = 0;
H1_int  = 0;        % mu neq 0
H0_int  = 0;        % mu = 0

% Pi Vector
l = linspace(0,2*pi,N);

disp('Running MC...')
for k = 1:M

    % Generate Noise
    pd = makedist('normal');
    pd.mu = 0;
    pd.sigma = 1;
    n0 = random(pd,N,1);
    
    % Inteference
    S = zeros(N,2);
    S(1:N,1)=cos(l.*(2*pi*f0));
    S(1:N,2)=sin(l.*(2*pi*f0));
    phi = A.*[cos(phi0);sin(phi0)];
    
    I = S*phi;
    
    % Create signal "917": Would add real signal model here
    theta = [(83/94);21/47;161/94];
    H = [1 -1 5;
        4 2 -2;
        7 -2 1];
    
    SIG = H*theta;
    
    % Observations
    if SIGNAL
       
            Y = SIG+I+n0;
    else
            Y = 0*SIG+I+n0;
        end
    
    % % Projections
    PHs = S*inv(S'*S)*S';
    PHc = eye(N)-PHs;
    PG = PHc*SIG*inv(SIG'*PHc*SIG)*SIG'*PHc;
    
    % Test statistics
%     Tkn = (Y'*PHs*Y)/((N)*(pd.sigma));
%     Tun = (Y'*PHs*Y)/(Y'*PHc*Y); %CFAR    
    Tint = (Y'*PG*Y)/(N*pd.sigma);
    
    if Tint>THRSH
        H1_int = H1_int+1;
    else
        H0_int = H0_int+1;
    end
%    
%     if Tkn>THRSH
%         H1_kn = H1_kn+1;
%     else
%         H0_kn = H0_kn+1;
%     end
% %     
%     if Tun>THRSH
%         H1_uk = H1_uk+1;
%     else
%         H0_uk = H0_uk+1;
%     end
end
% H1_kn
% H0_kn
% H1_uk
% H0_uk
H1_int
H0_int
