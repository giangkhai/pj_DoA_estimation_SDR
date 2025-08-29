import numpy as np
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.linalg import toeplitz

def resolve_phys_arr(d_Nr, d_P, d_SNR_low, d_SNR_high, d_ang, d_N, d_beta , d_num_trials): 
    # ==========================================================================================
    # --- QPSK generator ---
    def s_qpsk(d_P):
        v_symbols = np.array([np.exp(1j * np.pi / 4), np.exp(3j * np.pi / 4),
                              np.exp(5j * np.pi / 4), np.exp(-1j * np.pi / 4)])
        return np.random.choice(v_symbols, size=(1, d_P))
    # ==========================================================================================
    # --- Coupling matrix ---
    def m_coupling(d_Nr, d_N, d_beta=0.2):
        B = min(d_Nr, d_N)
        v_c = np.zeros(B + 1, dtype=complex)
        v_c[0] = 1.0
        for l in range(1, B + 1):
            v_c[l] = (d_beta / l) * np.exp(-2 * np.pi * l / (B + 1))
        v_coeffs = np.concatenate([v_c, np.zeros(d_Nr - (B + 1))])
        return toeplitz(v_coeffs)
    m_C = m_coupling(d_Nr, d_N, d_beta)
    # ==========================================================================================
    # --- Spectrum normalization ---
    def normalize_spectrum(v_spec):
        v_spec = np.array(v_spec) - np.max(v_spec)
        return 10 * np.log10(10 ** (v_spec / 10) / np.max(10 ** (v_spec / 10)))
    # ==========================================================================================
    # --- Get top angles ---
    def get_top_angles(v_peaks, v_spec, d_Ne, v_theta_g):
        if len(v_peaks) == 0: 
            return []
        v_peak_heights = v_spec[v_peaks]
        v_top_peaks = v_peaks[np.argsort(v_peak_heights)[::-1][:d_Ne]]
        return sorted([-90 + 180 * i / (len(v_theta_g) - 1) for i in v_top_peaks])
    # ==========================================================================================
    # --- Detection all-or-none ---
    def detection_ratio_all_or_none(v_true_angles, v_detected_angles, d_tol=1.0):
        v_true_angles = v_true_angles.flatten()
        all_detected = all(any(np.abs(angle - np.array(v_detected_angles)) <= d_tol) for angle in v_true_angles)
        return 1.0 if all_detected else 0.0
    # ==========================================================================================
    # --- Initialize ---
    v_theta_deg = np.array([d_ang]).reshape(1, -1)
    v_theta_rad = np.deg2rad(v_theta_deg)
    v_theta_g = np.linspace(-np.pi / 2, np.pi / 2, 720)
    v_SNR_range = np.arange(d_SNR_low, d_SNR_high + 1, 1)  
    v_ratios_MUSIC = []
    # ==========================================================================================
    # --- Main loop over SNRs ---
    for d_SNR in tqdm(v_SNR_range, desc="Processing SNRs"):
        d_sigma = 10 ** (-d_SNR / 20)
        d_success_MUSIC = 0

        for _ in range(d_num_trials):
            # Steering matrix and signals
            m_A = np.exp(-2j * np.pi * 0.5 * np.arange(d_Nr).reshape(-1, 1) @ np.sin(v_theta_rad).reshape(1, -1))
            m_s = np.vstack([s_qpsk(d_P) for _ in range(1)])  # single source
            m_y = m_C @ m_A @ m_s
            m_n = np.random.randn(d_Nr, d_P) + 1j * np.random.randn(d_Nr, d_P)
            m_y += d_sigma / np.sqrt(2) * m_n

            # Covariance
            m_R = np.cov(m_y, rowvar=True, bias=True)

            # Eigen decomposition
            v_eigvals, m_eigvecs = np.linalg.eig(m_R)
            v_eig_order = np.argsort(np.abs(v_eigvals))
            m_V = m_eigvecs[:, v_eig_order[:d_Nr - 1]]  # Noise subspace

            # MUSIC spectrum
            v_results_MUSIC = []
            for theta_i in v_theta_g:
                v_w = np.exp(-2j * np.pi * 0.5 * np.arange(d_Nr) * np.sin(theta_i)).reshape(-1, 1)
                d_P_MUSIC = (1 / (v_w.conj().T @ m_V @ m_V.conj().T @ v_w)).item()
                v_results_MUSIC.append(10 * np.log10(np.abs(d_P_MUSIC)))

            v_results_MUSIC = normalize_spectrum(v_results_MUSIC)
            v_top_angles = get_top_angles(find_peaks(v_results_MUSIC, height=-100)[0], v_results_MUSIC, 1, v_theta_g)
            d_success_MUSIC += detection_ratio_all_or_none(v_theta_deg, v_top_angles)

        v_ratios_MUSIC.append(d_success_MUSIC / d_num_trials)

    return v_SNR_range, v_ratios_MUSIC
# ============================================================================================================================
# ===== Example  =====
""" Only consider, investigate the case of single-source """     
d_Nr = 5              # Number of antennas
d_P = 100             # Number of snapshots 
d_SNR_low = -10       # lower bound of tested SNR
d_SNR_high = 0        # upper bound of tested SNR 
d_ang = 10            # angle 1 in degrees
d_N = 4               # Banded of coupling coefficient
d_beta = 0.01         # Largest agnitude of coefficient
d_num_trials = 10     # Number of Monte Carlo runs
v_SNR_range, v_ratios_MUSIC = resolve_phys_arr(d_Nr, d_P, d_SNR_low, d_SNR_high, d_ang, d_N, d_beta , d_num_trials)
# ============================================================================================================================
plt.figure(figsize=(8, 5))
plt.plot(v_SNR_range, v_ratios_MUSIC, '^-', label='MUSIC')
plt.xlabel('SNR (dB)')
plt.ylabel('Detection probability (all sources detected)')
plt.title(f'Detection probability vs SNR ({d_num_trials} trials per SNR)')
plt.grid(True)
plt.ylim(0, 1.05)
plt.legend()
plt.show()
