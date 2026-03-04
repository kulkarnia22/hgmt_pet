#!/usr/bin/env python3
"""
Upgraded list-mode PET FBP with:
  - "manual MSRB" for continuous-z detectors (virtual axial slices + weighted multi-slice deposition)
  - plane projection: per-slice sinogram uses intersection of the 3D LOR with z=z_k planes
  - bilinear deposition in BOTH phi and s (greatly reduces star/streak artifacts)
  - optional windowed ramp filter (debug), default is pure Ram–Lak
  - .voxels output compatible with your existing render.py (same header + one frame)

Input LOR file (binary):
  per event: 6 float64 values (x1 y1 z1 x2 y2 z2) in mm

Output .voxels:
  float64 X_LENGTH
  float64 Y_LENGTH
  float64 Z_LENGTH
  int32   iterations (we write 1)
  int32   X_RES
  int32   Y_RES
  int32   Z_RES
  then float64 volume values for 1 frame, stored in (X_RES, Y_RES, Z_RES) with Z fastest.

Notes / recommendations:
  - For PSF debugging, reconstruct a SMALL transverse ROI around the source (e.g. R_IMG=10 mm)
    then later expand to R_IMG=50 mm if you want.
  - This is still analytic FBP (not TOF, no attenuation, no normalization).
  - For NEMA-like "FBP + ramp only", set FILTER_WINDOW="none".

Run:
  python fbp_upgraded.py
View:
  python render.py image_fbp.voxels
"""

import numpy as np

# ===============================
# USER PARAMETERS (hard-coded)
# ===============================

# Input / output
LOR_FILE = "data/HGMTPointVac_no_tof.lor"   # 6 float64 per event: x1 y1 z1 x2 y2 z2 (mm)
OUT_FILE = "data/image_fbp.voxels"
DTYPE = np.float64                   # must match how you wrote the file (double -> float64)

# Reconstruction grid (mm)
# For PSF debugging, use R_IMG=5..10 mm. Later you can increase to 50 mm.
R_IMG = 10     # transverse half-FOV (mm). Reconstruct x,y in [-R_IMG, +R_IMG]
DX = 0.05        # pixel size (mm) in x/y

# Axial slab (mm) to reconstruct (you can keep small for point source)
Z_MIN = -10.0
Z_MAX =  10.0
DZ = 0.5

# Sinogram sampling
N_PHI = 360      # angles in [0, pi). 180 is a good start; 360 is ok once stable.
DS = DX          # s-bin spacing (mm)
S_MAX = R_IMG    # max |s| (mm) for sinogram

# "Manual MSRB" parameters (continuous-z)
MSRB_HALF_SLICES = 2          # each LOR contributes to k0-MS..k0+MS (so 2 => 5 slices)
MSRB_KERNEL = "tri"           # "tri" or "gauss"
MSRB_WIDTH_MM = 2.0 * DZ      # kernel width parameter in mm:
                              # tri: weights go to 0 at |dz|=MSRB_WIDTH_MM
                              # gauss: sigma = MSRB_WIDTH_MM

# Filter
FILTER_WINDOW = "none"        # "none" (pure ramp), "hann" (debug), "shepp" (debug)

# ===============================
# Internal / derived parameters
# ===============================

EPS = 1e-12

# Derived sinogram bins
N_S = int(np.ceil((2.0 * S_MAX) / DS))
if N_S % 2 == 1:
    N_S += 1
DS_EFF = (2.0 * S_MAX) / N_S
S0 = -S_MAX + 0.5 * DS_EFF   # first bin center

# Derived image grid
x_centers = np.arange(-R_IMG + DX / 2, R_IMG - DX / 2 + EPS, DX, dtype=np.float64)
y_centers = np.arange(-R_IMG + DX / 2, R_IMG - DX / 2 + EPS, DX, dtype=np.float64)
X_RES = x_centers.size
Y_RES = y_centers.size

z_edges = np.arange(Z_MIN, Z_MAX + DZ, DZ, dtype=np.float64)
z_centers = 0.5 * (z_edges[:-1] + z_edges[1:])
Z_RES = z_centers.size

phis = (np.arange(N_PHI, dtype=np.float64) + 0.5) * (np.pi / N_PHI)
DPHI = np.pi / N_PHI


# ===============================
# Helpers
# ===============================

def wrap_phi_0_pi(phi: np.ndarray) -> np.ndarray:
    return np.mod(phi, np.pi)


def ramp_filter_fft(proj: np.ndarray, ds: float, window: str = "none") -> np.ndarray:
    """
    Ramp (Ram–Lak) filter in frequency domain; optional apodization for debugging.
    proj: (N_S,)
    ds: spacing in s (mm)

    NOTE: Absolute scaling doesn't matter for PSF width; relative shape does.
    """
    n = proj.shape[0]
    freqs = np.fft.rfftfreq(n, d=ds)          # cycles/mm
    w = np.abs(freqs)                         # ramp

    if window == "hann":
        # Hann on [0, f_nyq]
        fmax = freqs[-1] if freqs.size else 0.0
        if fmax > 0:
            win = 0.5 * (1.0 + np.cos(np.pi * freqs / fmax))
            w = w * win
    elif window == "shepp":
        # Shepp-Logan window: sinc(pi f / (2 fmax))
        fmax = freqs[-1] if freqs.size else 0.0
        if fmax > 0:
            x = freqs / fmax
            # avoid division by zero at f=0: sinc(0)=1
            win = np.sinc(x / 2.0)
            w = w * win
    elif window == "none":
        pass
    else:
        raise ValueError(f"Unknown FILTER_WINDOW={window}")

    F = np.fft.rfft(proj)
    return np.fft.irfft(F * w, n=n)


def backproject_slice(sino_filt: np.ndarray,
                      cos_phi: np.ndarray,
                      sin_phi: np.ndarray,
                      x_centers: np.ndarray,
                      y_centers: np.ndarray,
                      s0: float,
                      ds: float) -> np.ndarray:
    """
    Backproject filtered sinogram into a 2D image.
    sino_filt: (N_PHI, N_S)
    Returns: (Y_RES, X_RES) for imshow
    """
    Xg, Yg = np.meshgrid(x_centers, y_centers, indexing="xy")  # (Y, X)
    img = np.zeros((y_centers.size, x_centers.size), dtype=np.float64)

    inv_ds = 1.0 / ds
    n_s = sino_filt.shape[1]

    for p in range(sino_filt.shape[0]):
        s = Xg * cos_phi[p] + Yg * sin_phi[p]
        u = (s - s0) * inv_ds
        i0 = np.floor(u).astype(np.int64)
        w = u - i0

        mask = (i0 >= 0) & (i0 < (n_s - 1))
        if not np.any(mask):
            continue

        proj = sino_filt[p]
        vals = np.zeros_like(img)
        vals[mask] = (1.0 - w[mask]) * proj[i0[mask]] + w[mask] * proj[i0[mask] + 1]
        img += vals * DPHI

    return img


def write_voxels(path: str, vol_zyx: np.ndarray, x_len: float, y_len: float, z_len: float) -> None:
    """
    Write .voxels file compatible with your render.py / EM output format.

    vol_zyx: shape (Z_RES, Y_RES, X_RES)
    stored as (X_RES, Y_RES, Z_RES) with Z fastest in file.
    """
    vol_zyx = np.asarray(vol_zyx, dtype=np.float64)
    z_res, y_res, x_res = vol_zyx.shape
    vol_xyz = np.transpose(vol_zyx, (2, 1, 0))  # (X, Y, Z)

    with open(path, "wb") as f:
        np.array([x_len], dtype=np.float64).tofile(f)
        np.array([y_len], dtype=np.float64).tofile(f)
        np.array([z_len], dtype=np.float64).tofile(f)
        np.array([np.int32(1)], dtype=np.int32).tofile(f)   # iterations = 1
        np.array([np.int32(x_res)], dtype=np.int32).tofile(f)
        np.array([np.int32(y_res)], dtype=np.int32).tofile(f)
        np.array([np.int32(z_res)], dtype=np.int32).tofile(f)
        vol_xyz.ravel(order="C").tofile(f)


def msrb_weights(dz: np.ndarray, kernel: str, width: float) -> np.ndarray:
    """
    Compute nonnegative weights as a function of dz = z_k - z_ref (mm).
    Returns unnormalized weights (caller should normalize per-event).
    """
    adz = np.abs(dz)
    if kernel == "tri":
        # 1 at dz=0, 0 at |dz|=width
        w = np.maximum(0.0, 1.0 - adz / max(width, EPS))
    elif kernel == "gauss":
        sigma = max(width, EPS)
        w = np.exp(-0.5 * (adz / sigma) ** 2)
    else:
        raise ValueError(f"Unknown MSRB_KERNEL={kernel}")
    return w


# ===============================
# Main
# ===============================

def main():
    # ---------- Read LORs ----------
    raw = np.fromfile(LOR_FILE, dtype=DTYPE)
    if raw.size % 6 != 0:
        raise ValueError(f"{LOR_FILE}: {raw.size} numbers not divisible by 6.")
    lor = raw.reshape((-1, 6))
    x1, y1, z1, x2, y2, z2 = lor.T
    n_total = lor.shape[0]
    print(f"Loaded {n_total} LORs from {LOR_FILE}")

    # ---------- Precompute transverse direction + normal (phi) ----------
    dx12 = x2 - x1
    dy12 = y2 - y1
    Lxy = np.sqrt(dx12 * dx12 + dy12 * dy12)
    good = Lxy > EPS

    x1 = x1[good]; y1 = y1[good]; z1 = z1[good]
    x2 = x2[good]; y2 = y2[good]; z2 = z2[good]
    dx12 = dx12[good]; dy12 = dy12[good]; Lxy = Lxy[good]
    n_use = x1.size
    print(f"Using {n_use} LORs after dropping near-zero XY-length events.")

    # normal in XY
    nx = -dy12 / Lxy
    ny =  dx12 / Lxy

    phi = wrap_phi_0_pi(np.arctan2(ny, nx))  # [0, pi)
    # continuous phi-bin coordinate in [0, N_PHI)
    u_phi = (phi / np.pi) * N_PHI
    i_phi0 = np.floor(u_phi).astype(np.int64)
    w_phi = u_phi - i_phi0
    i_phi0 = np.mod(i_phi0, N_PHI)
    i_phi1 = np.mod(i_phi0 + 1, N_PHI)

    # precompute cos/sin for backprojection
    cos_phi = np.cos(phis)
    sin_phi = np.sin(phis)

    # ---------- Setup volume ----------
    vol_zyx = np.zeros((Z_RES, Y_RES, X_RES), dtype=np.float64)

    # ---------- Reconstruct slice-by-slice ----------
    print(f"Reconstruction grid: X_RES={X_RES}, Y_RES={Y_RES}, Z_RES={Z_RES} | DX={DX} mm, DZ={DZ} mm | ROI R={R_IMG} mm")
    print(f"Sinogram grid: N_PHI={N_PHI}, N_S={N_S}, ds={DS_EFF:.6f} mm, s in [-{S_MAX}, +{S_MAX}] mm")
    print(f"MSRB: half_slices={MSRB_HALF_SLICES} (=> {2*MSRB_HALF_SLICES+1} slices/LOR), kernel={MSRB_KERNEL}, width={MSRB_WIDTH_MM} mm")
    print(f"Filter window: {FILTER_WINDOW}")

    # We'll accumulate sinograms per slice in a streaming way:
    # For each LOR, choose a small neighborhood of slices around z_ref (midpoint),
    # then plane-project intersection into those slice planes, and deposit into those sinograms.
    # This avoids an N_LOR * Z_RES loop.
    sinos = np.zeros((Z_RES, N_PHI, N_S), dtype=np.float64)

    z_ref = 0.5 * (z1 + z2)

    # nearest slice index to z_ref
    k0 = np.rint((z_ref - z_centers[0]) / DZ).astype(np.int64)
    k0 = np.clip(k0, 0, Z_RES - 1)

    # For each offset in [-MS..MS], process those contributions
    for dk in range(-MSRB_HALF_SLICES, MSRB_HALF_SLICES + 1):
        kk = k0 + dk
        valid_k = (kk >= 0) & (kk < Z_RES)
        if not np.any(valid_k):
            continue

        kk_v = kk[valid_k]
        z1_v = z1[valid_k]; z2_v = z2[valid_k]
        x1_v = x1[valid_k]; y1_v = y1[valid_k]
        x2_v = x2[valid_k]; y2_v = y2[valid_k]
        nx_v = nx[valid_k]; ny_v = ny[valid_k]
        i0_v = i_phi0[valid_k]; i1_v = i_phi1[valid_k]; wphi_v = w_phi[valid_k]
        zref_v = z_ref[valid_k]

        z_k = z_centers[kk_v]  # slice plane z values per event

        # plane intersection parameter t = (z_k - z1)/(z2 - z1)
        denom = (z2_v - z1_v)
        good_denom = np.abs(denom) > EPS
        if not np.any(good_denom):
            continue

        kk_v = kk_v[good_denom]
        z_k = z_k[good_denom]
        x1_v = x1_v[good_denom]; y1_v = y1_v[good_denom]; z1_v = z1_v[good_denom]
        x2_v = x2_v[good_denom]; y2_v = y2_v[good_denom]; z2_v = z2_v[good_denom]
        nx_v = nx_v[good_denom]; ny_v = ny_v[good_denom]
        i0_v = i0_v[good_denom]; i1_v = i1_v[good_denom]; wphi_v = wphi_v[good_denom]
        zref_v = zref_v[good_denom]
        denom = denom[good_denom]

        t = (z_k - z1_v) / denom
        # Only if intersection lies between endpoints
        hit = (t >= 0.0) & (t <= 1.0)
        if not np.any(hit):
            continue

        kk_v = kk_v[hit]
        z_k = z_k[hit]
        x1_v = x1_v[hit]; y1_v = y1_v[hit]
        x2_v = x2_v[hit]; y2_v = y2_v[hit]
        nx_v = nx_v[hit]; ny_v = ny_v[hit]
        i0_v = i0_v[hit]; i1_v = i1_v[hit]; wphi_v = wphi_v[hit]
        zref_v = zref_v[hit]
        t = t[hit]

        # Intersection point (x_k, y_k)
        xk = x1_v + t * (x2_v - x1_v)
        yk = y1_v + t * (y2_v - y1_v)

        # Sinogram coordinate for this slice: s_k = n · (xk, yk)
        s = nx_v * xk + ny_v * yk

        # Drop events outside sinogram s range
        u_s = (s - S0) / DS_EFF
        is0 = np.floor(u_s).astype(np.int64)
        ws = u_s - is0
        ok_s = (is0 >= 0) & (is0 < (N_S - 1))
        if not np.any(ok_s):
            continue

        kk_v = kk_v[ok_s]
        i0_v = i0_v[ok_s]; i1_v = i1_v[ok_s]; wphi_v = wphi_v[ok_s]
        is0 = is0[ok_s]; ws = ws[ok_s]
        zref_v = zref_v[ok_s]
        z_k = z_k[ok_s]

        # MSRB weight based on dz from reference
        dz = (z_k - zref_v)
        wk = msrb_weights(dz, MSRB_KERNEL, MSRB_WIDTH_MM)

        # We'll normalize weights across dk later by accumulating a per-event norm.
        # But since we're processing dk separately and we dropped some events by intersection/s-range,
        # true per-event normalization is messy in vector form.
        #
        # Pragmatic approach (works well for point-source PSF):
        # use unnormalized wk here, then at the end normalize each slice by total deposited weight.
        #
        # Alternative: do a 2-pass normalization per event. That's heavier.
        #
        # We'll do slice-wise normalization via a sensitivity-like map in sinogram space (weight sums).
        # This keeps relative geometry correct and avoids bias from varying # of intersected slices.

        # Bilinear deposition into (phi, s):
        w00 = (1.0 - wphi_v) * (1.0 - ws) * wk
        w01 = (1.0 - wphi_v) * (ws)       * wk
        w10 = (wphi_v)       * (1.0 - ws) * wk
        w11 = (wphi_v)       * (ws)       * wk

        # Deposit into sinogram of slice kk_v
        # phi bin i0_v/i1_v, s bin is0/is0+1
        np.add.at(sinos, (kk_v, i0_v, is0),     w00)
        np.add.at(sinos, (kk_v, i0_v, is0 + 1), w01)
        np.add.at(sinos, (kk_v, i1_v, is0),     w10)
        np.add.at(sinos, (kk_v, i1_v, is0 + 1), w11)

    # ---------- Filter + backproject each slice ----------
    for kk in range(Z_RES):
        sino = sinos[kk]

        # If slice has essentially no counts, skip
        total_counts = np.sum(sino)
        if total_counts <= 0:
            continue

        # Normalize sinogram by its total "weight" so slices are comparable in scale
        # (doesn't affect FWHM, only intensity scaling)
        sino = sino * (1.0 / max(total_counts, EPS))

        # Filter each projection
        sino_filt = np.empty_like(sino)
        for p in range(N_PHI):
            sino_filt[p] = ramp_filter_fft(sino[p], DS_EFF, window=FILTER_WINDOW)

        # Backproject
        img_yx = backproject_slice(sino_filt, cos_phi, sin_phi, x_centers, y_centers, S0, DS_EFF)
        vol_zyx[kk] = img_yx

        if kk % max(1, Z_RES // 10) == 0:
            print(f"Reconstructed slice {kk+1}/{Z_RES} (z={z_centers[kk]:.3f} mm), sinogram weight sum={total_counts:.3e}")

    # ---------- Write output ----------
    x_len = 2.0 * R_IMG
    y_len = 2.0 * R_IMG
    z_len = (Z_MAX - Z_MIN)
    write_voxels(OUT_FILE, vol_zyx, x_len=x_len, y_len=y_len, z_len=z_len)

    print(f"\nSaved upgraded FBP volume to: {OUT_FILE}")
    print(f"View with: python render.py {OUT_FILE}")
    print("\nIf the image is still streaky:")
    print("  - Reduce N_PHI to 90 or 120 for debugging")
    print("  - Keep R_IMG small (5-10 mm) until geometry looks right")
    print("  - Try FILTER_WINDOW='hann' to verify the PSF becomes compact (debug only)")
    print("  - Then set FILTER_WINDOW back to 'none' for ramp-only comparisons")


if __name__ == "__main__":
    main()