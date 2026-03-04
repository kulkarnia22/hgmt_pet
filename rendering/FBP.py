#!/usr/bin/env python3
"""
List-mode 2D FBP (slice-by-slice) for PET coincidences, with .voxels output
compatible with your render.py / imager.c format.

Input LOR file format (binary):
  Per event: 6 float64 (double) values in mm:
    x1 y1 z1 x2 y2 z2

Output .voxels format (binary), matching your EM output:
  float64 X_LENGTH
  float64 Y_LENGTH
  float64 Z_LENGTH
  int32   iterations   (for FBP we write 1 frame)
  int32   X_RES
  int32   Y_RES
  int32   Z_RES
  then float64 image values of length iterations * X_RES * Y_RES * Z_RES
  stored in (x, y, z) order with z fastest, i.e. C-order for array shape (X_RES, Y_RES, Z_RES).

Run:
  python fbp.py
Then view:
  python render.py image_fbp.voxels
"""

import numpy as np

# ===============================
# Hard-coded reconstruction setup
# ===============================

# Input / output
LOR_FILE = "data/HGMTPointVac_no_tof.lor"     # binary: 6 float64 per event: x1 y1 z1 x2 y2 z2 (mm)
OUT_FILE = "data/image_point_vac_fbp.voxels"         # output compatible with render.py
DTYPE = np.float64                    # matches C 'double'

# Reconstruction grid (mm)
R_IMG = 50.0      # transverse half-FOV radius; reconstruct x,y in [-R_IMG, +R_IMG]
DX = 0.05         # pixel size in x/y
Z_MIN = -5.0      # reconstruction slab min z
Z_MAX = 5.0       # reconstruction slab max z
DZ = 0.5          # slice thickness

# Tilt gating (minimize oblique LORs for 2D approximation)
#Z_TILT_MAX = 2.0 * DZ   # keep if |z2 - z1| <= Z_TILT_MAX; try DZ, 2*DZ, 4*DZ, np.inf
Z_TILT_MAX = np.inf

# Sinogram sampling
N_PHI = 360       # number of angles in [0, pi)
S_MAX = R_IMG     # max |s| in mm
DS = DX           # s-bin spacing (mm) ~ pixel size is a good start
N_S = int(np.ceil((2.0 * S_MAX) / DS))
if N_S % 2 == 1:
    N_S += 1      # make even for FFT convenience

# Linear deposition in sinogram s direction (reduces binning blur)
LINEAR_DEPOSIT = True  # set False to use nearest-bin only

# -----------------------------
# Internal constants
# -----------------------------
EPS = 1e-12


def wrap_phi_0_pi(phi: np.ndarray) -> np.ndarray:
    """Wrap angles to [0, pi)."""
    return np.mod(phi, np.pi)


def ramp_filter_fft(proj: np.ndarray, ds: float) -> np.ndarray:
    """
    Ram–Lak (ramp) filter with no apodization/windowing.
    proj: shape (N_S,)
    ds: spacing in s (mm)

    Overall amplitude scaling is not important for PSF/FWHM comparisons.
    """
    n = proj.shape[0]
    freqs = np.fft.rfftfreq(n, d=ds)  # cycles/mm
    ramp = np.abs(freqs)
    F = np.fft.rfft(proj)
    return np.fft.irfft(F * ramp, n=n)


def backproject_slice(
    sino_filt: np.ndarray,
    phis: np.ndarray,
    s0: float,
    ds: float,
    x_centers: np.ndarray,
    y_centers: np.ndarray,
) -> np.ndarray:
    """
    Backproject a filtered sinogram into a 2D image.

    sino_filt: (N_PHI, N_S)
    phis:      (N_PHI,)
    s0:        first s-bin center
    ds:        s-bin spacing
    x_centers: (X_RES,)
    y_centers: (Y_RES,)

    Returns img with shape (Y_RES, X_RES) for imshow.
    """
    Xg, Yg = np.meshgrid(x_centers, y_centers, indexing="xy")  # (Y_RES, X_RES)
    img = np.zeros((y_centers.size, x_centers.size), dtype=np.float64)

    dphi = np.pi / phis.size
    inv_ds = 1.0 / ds
    cos_phi = np.cos(phis)
    sin_phi = np.sin(phis)

    n_s = sino_filt.shape[1]

    for i in range(phis.size):
        s = Xg * cos_phi[i] + Yg * sin_phi[i]
        u = (s - s0) * inv_ds
        i0 = np.floor(u).astype(np.int64)
        w = u - i0

        mask = (i0 >= 0) & (i0 < (n_s - 1))
        if not np.any(mask):
            continue

        proj = sino_filt[i]
        vals = np.zeros_like(img)
        vals[mask] = (1.0 - w[mask]) * proj[i0[mask]] + w[mask] * proj[i0[mask] + 1]

        img += vals * dphi

    return img


def write_voxels(path: str, vol_zyx: np.ndarray, x_len: float, y_len: float, z_len: float) -> None:
    """
    Write .voxels file compatible with your render.py and imager.c output format.

    vol_zyx: shape (Z_RES, Y_RES, X_RES), float64
    Writes 1 frame, stored as (X_RES, Y_RES, Z_RES) in C order.
    """
    vol_zyx = np.asarray(vol_zyx, dtype=np.float64)
    if vol_zyx.ndim != 3:
        raise ValueError(f"Expected 3D volume, got shape {vol_zyx.shape}")

    z_res, y_res, x_res = vol_zyx.shape
    vol_xyz = np.transpose(vol_zyx, (2, 1, 0))  # (X_RES, Y_RES, Z_RES)

    iterations = np.int32(1)
    with open(path, "wb") as f:
        np.array([x_len], dtype=np.float64).tofile(f)
        np.array([y_len], dtype=np.float64).tofile(f)
        np.array([z_len], dtype=np.float64).tofile(f)
        np.array([iterations], dtype=np.int32).tofile(f)
        np.array([x_res], dtype=np.int32).tofile(f)
        np.array([y_res], dtype=np.int32).tofile(f)
        np.array([z_res], dtype=np.int32).tofile(f)
        vol_xyz.ravel(order="C").tofile(f)


def main():
    # -----------------------------
    # Read LORs
    # -----------------------------
    data = np.fromfile(LOR_FILE, dtype=DTYPE)
    if data.size % 6 != 0:
        raise ValueError(f"{LOR_FILE}: file has {data.size} numbers, not divisible by 6 (expected 6 per LOR).")
    lor = data.reshape((-1, 6))
    x1, y1, z1, x2, y2, z2 = lor.T
    n_total = lor.shape[0]

    # -----------------------------
    # Grid setup
    # -----------------------------
    # x/y pixel centers (mm)
    x_centers = np.arange(-R_IMG + DX / 2, R_IMG - DX / 2 + EPS, DX, dtype=np.float64)
    y_centers = np.arange(-R_IMG + DX / 2, R_IMG - DX / 2 + EPS, DX, dtype=np.float64)
    x_res = x_centers.size
    y_res = y_centers.size

    # z slice centers (mm)
    z_edges = np.arange(Z_MIN, Z_MAX + DZ, DZ, dtype=np.float64)
    z_centers = 0.5 * (z_edges[:-1] + z_edges[1:])
    z_res = z_centers.size

    # sinogram bins
    phis = (np.arange(N_PHI, dtype=np.float64) + 0.5) * (np.pi / N_PHI)
    s_bins = (np.arange(N_S, dtype=np.float64) + 0.5) * (2.0 * S_MAX / N_S) - S_MAX
    s0 = s_bins[0]
    ds_eff = s_bins[1] - s_bins[0]

    # -----------------------------
    # Slice assignment + tilt gate (SSRB-like)
    # -----------------------------
    z_mid = 0.5 * (z1 + z2)
    dz_lor = np.abs(z2 - z1)

    k = np.floor((z_mid - Z_MIN) / DZ).astype(np.int64)
    in_range = (k >= 0) & (k < z_res)
    tilt_ok = dz_lor <= Z_TILT_MAX
    keep = in_range & tilt_ok

    x1 = x1[keep]; y1 = y1[keep]
    x2 = x2[keep]; y2 = y2[keep]
    k = k[keep]
    n_kept = k.size

    print(f"Loaded {n_total} LORs; kept {n_kept} after slice-range + tilt gate (|z2-z1| <= {Z_TILT_MAX:.3f} mm).")
    print(f"Reconstruction grid: X_RES={x_res}, Y_RES={y_res}, Z_RES={z_res}  | DX={DX} mm, DZ={DZ} mm")
    print(f"Sinogram grid: N_PHI={N_PHI}, N_S={N_S}, ds={ds_eff:.6f} mm, s in [-{S_MAX}, +{S_MAX}] mm")

    # -----------------------------
    # Compute (s, phi) from transverse endpoints
    # -----------------------------
    dx12 = x2 - x1
    dy12 = y2 - y1
    L = np.sqrt(dx12 * dx12 + dy12 * dy12)
    good_len = L > EPS

    x1 = x1[good_len]; y1 = y1[good_len]
    dx12 = dx12[good_len]; dy12 = dy12[good_len]
    k = k[good_len]
    L = L[good_len]

    # unit normal n = (-dy, dx)/L
    nx = -dy12 / L
    ny =  dx12 / L

    s = nx * x1 + ny * y1
    phi = wrap_phi_0_pi(np.arctan2(ny, nx))

    # phi bins
    i_phi = np.floor(phi / np.pi * N_PHI).astype(np.int64)
    i_phi = np.clip(i_phi, 0, N_PHI - 1)

    # s bins
    u = (s - s0) / ds_eff
    if LINEAR_DEPOSIT:
        i_s0 = np.floor(u).astype(np.int64)
        w = u - i_s0
        valid_s = (i_s0 >= 0) & (i_s0 < (N_S - 1))
        k = k[valid_s]; i_phi = i_phi[valid_s]; i_s0 = i_s0[valid_s]; w = w[valid_s]
    else:
        i_s0 = np.rint(u).astype(np.int64)
        valid_s = (i_s0 >= 0) & (i_s0 < N_S)
        k = k[valid_s]; i_phi = i_phi[valid_s]; i_s0 = i_s0[valid_s]
        w = None

    # -----------------------------
    # Allocate volume (Z, Y, X)
    # -----------------------------
    vol_zyx = np.zeros((z_res, y_res, x_res), dtype=np.float64)

    # -----------------------------
    # Reconstruct slice-by-slice
    # -----------------------------
    for kk in range(z_res):
        idx = np.where(k == kk)[0]
        if idx.size == 0:
            continue

        # Build sinogram counts (phi, s)
        sino = np.zeros((N_PHI, N_S), dtype=np.float64)
        ip = i_phi[idx]
        is0 = i_s0[idx]

        if LINEAR_DEPOSIT:
            ww = w[idx]
            np.add.at(sino, (ip, is0),     (1.0 - ww))
            np.add.at(sino, (ip, is0 + 1), (ww))
        else:
            np.add.at(sino, (ip, is0), 1.0)

        # Filter each projection with pure Ram–Lak
        sino_filt = np.empty_like(sino)
        for p in range(N_PHI):
            sino_filt[p] = ramp_filter_fft(sino[p], ds_eff)

        # Backproject
        img_yx = backproject_slice(sino_filt, phis, s0, ds_eff, x_centers, y_centers)
        vol_zyx[kk] = img_yx

        if kk % max(1, z_res // 10) == 0:
            print(f"Reconstructed slice {kk+1}/{z_res} (z={z_centers[kk]:.3f} mm) with {idx.size} events.")

    # -----------------------------
    # Write .voxels output
    # -----------------------------
    x_len = 2.0 * R_IMG
    y_len = 2.0 * R_IMG
    z_len = (Z_MAX - Z_MIN)

    write_voxels(OUT_FILE, vol_zyx, x_len=x_len, y_len=y_len, z_len=z_len)

    print(f"\nSaved FBP volume to: {OUT_FILE}")
    print(f"View with: python render.py {OUT_FILE}")


if __name__ == "__main__":
    main()