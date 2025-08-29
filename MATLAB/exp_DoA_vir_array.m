%% ===================================================================================================================================
%% This code is used for the case where antenna 1 is fixed and antenna 2 is moving.
%% When d_l = 1, it corresponds to using MUSIC for 2 physical antennas.
%% The recorded files must be named rx1_l and rx2_l, and they must have the same length for each acquisition, with l = 1, 2, ..., L.
%% ===================================================================================================================================
%% Parameters
d_n = 4;                          % d_n x d_n virtual array
d_a = 0.336;                      % antenna spacing
d_fc = 446.00625e6;               % center frequency
d_c = 299792458;
d_lambda = d_c/d_fc;
d_calib_phase = 0.985202754919110 + 0.171392916130733i; % phase calib
d_fs = 1e6;
%% ===================================================================================================================================
%% Read and form IQ signals
v_rx = cell(d_n-1,2); % cell{L,1} = rx1_L, cell{L,2} = rx2_L
for L = 1:d_n-1
    f1 = fopen(sprintf('rx1_%d.dat', L), 'rb');
    f2 = fopen(sprintf('rx2_%d.dat', L), 'rb');
    raw1 = fread(f1, 'float32'); fclose(f1);
    raw2 = fread(f2, 'float32'); fclose(f2);
    v_rx{L,1} = (raw1(1:2:end) + 1j*raw1(2:2:end)) * d_calib_phase;
    v_rx{L,2} = (raw2(1:2:end) + 1j*raw2(2:2:end));
end
%% ===================================================================================================================================
%% Compute cross-correlation and powers
m_R = zeros(d_n,d_n);
for i = 1:d_n
    for j = 1:d_n
        if i==j
            % diagonal: power of corresponding acquisition
            if i==1
                m_R(i,j) = mean(abs(v_rx{1,1}).^2);
            else
                m_R(i,j) = mean(abs(v_rx{i-1,2}).^2);
            end
        else
            L = max(i,j)-1;
            m_R(i,j) = mean(v_rx{L,1} .* conj(v_rx{L,2}));
        end
    end
end
%% ===================================================================================================================================
%% Steering vector function
v_w = @(theta_deg) exp(-1i*2*pi*d_a/d_lambda*(0:d_n-1).' * sind(theta_deg));
%% ===================================================================================================================================
%% MUSIC
[m_V, ~] = eig(m_R);
[~, idx] = sort(diag(m_R),'descend'); % or eigvals if using eig(R)
m_V = m_V(:, idx);
m_En = m_V(:, d_n:end); % Noise subspace
v_angles = linspace(-90,90,360);
v_P = zeros(size(v_angles));
for k = 1:length(v_angles)
    v_a = v_w(v_angles(k));
    v_P(k) = 1 / abs(v_a'*(m_En*m_En')*v_a);
end
%% ===================================================================================================================================
%% Peak detection
[v_pks, v_locs] = findpeaks(v_P,'SortStr','descend','NPeaks',d_n-1,'MinPeakProminence',1e-9);
v_est_doa = v_angles(sort(v_locs));
fprintf("Estimated DOAs: "); fprintf('%.4f ', v_est_doa); fprintf('\n');
