import numpy as np
from tqdm import tqdm
from scipy.signal import find_peaks
import matplotlib.pyplot as plt

def vir_UCA(d_Nr: int, d_Nt: int, d_P: int, d_SNR: float, d_N_test: int = 1):
    # ================================================================================================
    # ===== Constants =====
    c_c = 3e8                                         # Speed of light
    c_fc = 39e6                                       # Carrier frequency 
    Ne = d_Nt                                         # Number of expected sources
    d_d = c_c / (8 * c_fc * np.sin(np.pi/d_Nr)**2)    # UCA antenna spacing
    # ================================================================================================
    # True DoA
    v_theta_deg_true = np.array([-60, -45, -30, 0, 15, 30, 45, 60])
    v_theta_rad_true = np.deg2rad(v_theta_deg_true)
    # ================================================================================================
    # QPSK generator
    def s_qpsk(d_P: int):
        v_const = np.array([np.exp(1j*np.pi/4), np.exp(3j*np.pi/4),
                            np.exp(5j*np.pi/4), np.exp(7j*np.pi/4)])
        return np.random.choice(v_const, size=(d_P,1))
    # ================================================================================================
    d_sigma = 10**(-d_SNR/10)
    v_err_list = []

    v_last_spectrum = None
    v_last_theta_deg = None

    for _ in tqdm(range(d_N_test), desc=f"Processing MUSIC SNR={d_SNR}dB", ncols=100, unit="iteration"):
        # Covariance matrix
        m_R = np.zeros((d_Nr,d_Nr), dtype=complex)
        for m_idx in range(d_Nr):
            for n_idx in range(d_Nr):
                s_m, s_mn = np.zeros((d_P,1), complex), np.zeros((d_P,1), complex)
                for theta in v_theta_rad_true:
                    s_base = s_qpsk(d_P)
                    phi_m  = 2*np.pi*c_fc*d_d/c_c * np.cos(theta - 2*np.pi*m_idx/d_Nr)
                    phi_mn = 2*np.pi*c_fc*d_d/c_c * np.cos(theta - 2*np.pi*(m_idx+n_idx)/d_Nr)
                    s_m  += s_base*np.exp(1j*phi_m)
                    s_mn += s_base*np.exp(1j*phi_mn)

                # Add noise
                s_m  += d_sigma/np.sqrt(2)*(np.random.randn(d_P,1)+1j*np.random.randn(d_P,1))
                s_mn += d_sigma/np.sqrt(2)*(np.random.randn(d_P,1)+1j*np.random.randn(d_P,1))

                if m_idx+n_idx < d_Nr and n_idx!=0:
                    m_R[m_idx,m_idx+n_idx] = np.mean(s_m*np.conj(s_mn))
                    m_R[m_idx+n_idx,m_idx] = np.conj(m_R[m_idx,m_idx+n_idx])
                elif n_idx==0:
                    m_R[m_idx,m_idx] = np.mean(s_m*np.conj(s_m))
        # ================================================================================================
        # MUSIC
        v_eigval, m_eigvec = np.linalg.eig(m_R)
        v_idx_order = np.argsort(np.abs(v_eigval))
        m_noise_space = m_eigvec[:, v_idx_order[:d_Nr-Ne]]

        v_theta_grid = np.linspace(-np.pi/2, np.pi/2, 720)
        v_spectrum = []
        for th in v_theta_grid:
            v_a = np.exp(2j*np.pi*c_fc*d_d/c_c * np.cos(th - 2*np.pi*np.arange(d_Nr)/d_Nr)).reshape(-1,1)
            d_Pm = 1 / (v_a.conj().T @ m_noise_space @ m_noise_space.conj().T @ v_a).item()
            v_spectrum.append(10*np.log10(np.abs(d_Pm)))
        v_spectrum = np.array(v_spectrum)
        # ================================================================================================
        # Peak detection
        v_peaks, _ = find_peaks(v_spectrum, height=5, distance=10, prominence=2)
        v_est_deg = [-90 + (180*i/(len(v_theta_grid)-1)) for i in v_peaks]
        v_est_deg = sorted(v_est_deg[:Ne])

        if len(v_est_deg) >= Ne:
            d_se = np.sum((np.array(v_est_deg)-v_theta_deg_true)**2)/Ne
            v_err_list.append(d_se)

        v_last_spectrum = v_spectrum - np.max(np.abs(v_spectrum))
        v_last_theta_deg = np.rad2deg(v_theta_grid)

    d_rmse = np.sqrt(np.mean(v_err_list)) if v_err_list else np.nan
    print(f"RMSE MUSIC @ {d_SNR} dB over {d_N_test} runs: {d_rmse:.4f}")

    return v_last_theta_deg, v_last_spectrum, d_rmse
# =========================================================================================================
# ===== Example =====
d_Nr = 12         # Number of virtual antennas
d_Nt = 8          # Number of sources
d_P = 2000        # Number of snapshots
d_SNR = 20        # Signal-to-noise ratio
d_N_test = 10     # Number of Monte Carlo runs

v_theta_grid, v_music_spec, d_rmse = vir_UCA(d_Nr, d_Nt, d_P, d_SNR, d_N_test)
# =========================================================================================================
# Plot last MUSIC spectrum
plt.figure(figsize=(10,6))
plt.plot(v_theta_grid, v_music_spec, 'b-', label="Normalized D-MUSIC Spectrum")
for ang in [-60, -45, -30, 0, 15, 30, 45, 60]:
    plt.axvline(x=ang, color='r', ls='--', alpha=0.5)
plt.title(f"MUSIC Spectrum @ {d_SNR} dB (last run)")
plt.xlabel("Angle (deg)")
plt.ylabel("Power (dB)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
