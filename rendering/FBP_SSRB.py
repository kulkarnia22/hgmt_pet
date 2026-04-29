"""
FBP Reconstruction for LOR endpoint data.
Produces a .voxels file readable by render.py.
 
Usage:
    python fbp_reconstruction.py <lor_endpoint_file> <output.voxels>
 
LOR file format: flat binary, 6 x float64 per LOR, no header.
Each record: x1 y1 z1 x2 y2 z2 (cm)
 
Output .voxels format:
    x_len  float64 (cm)
    y_len  float64 (cm)
    z_len  float64 (cm)
    iterations int32  (always 1 for FBP)
    x_res  int32
    y_res  int32
    z_res  int32
    x_res * y_res * z_res float64 voxel values (x-major order)
"""
 
import sys
import numpy as np
 
# ============================================================
# Configuration
# ============================================================
 
# Reconstruction volume half-extents in cm, centered on origin
# ±5 cm gives enough room to see the full PSF and surrounding background
RECON_X_CM = 5.0
RECON_Y_CM = 5.0
RECON_Z_CM = 5.0
 
# Voxel size in cm
# 0.05 cm = 0.5 mm — satisfies NEMA 1/3 FWHM requirement for 1-2 mm PSF
VOXEL_SIZE_CM = 0.05
 
# Sinogram parameters
N_ANGLES       = 360    # projection angles over [0, pi) — more angles = fewer streak artifacts
N_RADIAL       = 8192    # radial bins
RADIAL_HALF_CM = 5 # half-width of radial FOV in cm, covers detectors at 45-102.54 cm
 
# Only keep LORs whose axial tilt is below this angle
MAX_AXIAL_ANGLE_DEG = 10.0
 
# ============================================================
# Derived parameters — do not edit below this line
# ============================================================
 
MAX_AXIAL_TAN   = np.tan(np.radians(MAX_AXIAL_ANGLE_DEG))
 
X_LEN = 2.0 * RECON_X_CM
Y_LEN = 2.0 * RECON_Y_CM
Z_LEN = 2.0 * RECON_Z_CM
 
X_RES = int(round(X_LEN / VOXEL_SIZE_CM))
Y_RES = int(round(Y_LEN / VOXEL_SIZE_CM))
Z_RES = int(round(Z_LEN / VOXEL_SIZE_CM))
 
angles          = np.linspace(0, np.pi, N_ANGLES, endpoint=False)
r_edges         = np.linspace(-RADIAL_HALF_CM, RADIAL_HALF_CM, N_RADIAL + 1)
dr              = r_edges[1] - r_edges[0]
z_edges         = np.linspace(-RECON_Z_CM, RECON_Z_CM, Z_RES + 1)
slice_thickness = Z_LEN / Z_RES
 
# Pixel coordinate grids for backprojection — shape (X_RES, Y_RES)
px = np.linspace(-RECON_X_CM + VOXEL_SIZE_CM / 2,
                  RECON_X_CM - VOXEL_SIZE_CM / 2, X_RES)
py = np.linspace(-RECON_Y_CM + VOXEL_SIZE_CM / 2,
                  RECON_Y_CM - VOXEL_SIZE_CM / 2, Y_RES)
PX, PY = np.meshgrid(px, py, indexing='ij')  # (X_RES, Y_RES)
 
# Precompute cos and sin for all angles
COS_PHI = np.cos(angles)
SIN_PHI = np.sin(angles)
 
 
# ============================================================
# Step 1: Load LORs
# ============================================================
def read_lors(filepath):
    raw = np.fromfile(filepath, dtype=np.float64)
    if len(raw) % 6 != 0:
        raise ValueError(
            f"File length {len(raw)} not divisible by 6. "
            "Check LOR endpoint file format.")
    lors = raw.reshape(-1, 6)
    print(f"Loaded {len(lors):,} LORs")
    return lors
 
 
# ============================================================
# Step 2: Axial obliqueness cut
# ============================================================
def axial_cut(lors):
    dx = lors[:, 3] - lors[:, 0]
    dy = lors[:, 4] - lors[:, 1]
    dz = lors[:, 5] - lors[:, 2]
    transaxial = np.sqrt(dx**2 + dy**2)
    valid      = transaxial > 1e-9
    axial_tan  = np.zeros(len(lors))
    axial_tan[valid] = np.abs(dz[valid]) / transaxial[valid]
    mask = axial_tan < MAX_AXIAL_TAN
    print(f"Axial cut ({MAX_AXIAL_ANGLE_DEG:.1f} deg): "
          f"{mask.sum():,} / {len(lors):,} accepted "
          f"({100 * mask.mean():.1f}%)")
    return lors[mask]
 
 
# ============================================================
# Step 3: Bin into sinograms (SSRB) — consistent angle convention
# ============================================================
def build_sinograms(lors):
    x1, y1, z1 = lors[:, 0], lors[:, 1], lors[:, 2]
    x2, y2, z2 = lors[:, 3], lors[:, 4], lors[:, 5]
 
    dx         = x2 - x1
    dy         = y2 - y1
    transaxial = np.sqrt(dx**2 + dy**2)
    valid      = transaxial > 1e-9
 
    # Normalized LOR direction vector in transaxial plane
    lx = np.where(valid, dx / transaxial, 0.0)
    ly = np.where(valid, dy / transaxial, 0.0)
 
    # Signed perpendicular distance from origin to LOR
    # r = x1*ly - y1*lx
    # This is consistent with r = x*cos(phi) + y*sin(phi)
    # when phi = arctan2(lx, -ly) is the angle of the LOR normal (-ly, lx)
    r = np.where(valid, x1 * ly - y1 * lx, 0.0)
 
    # Angle of the LOR normal vector (-ly, lx), mapped to [0, pi)
    phi = np.arctan2(lx, -ly)
    neg = phi < 0
    phi = np.where(neg, phi + np.pi, phi)
    r   = np.where(neg, -r, r)
 
    # Axial midpoint for SSRB slice assignment
    z_mid = 0.5 * (z1 + z2)
 
    # Bin indices
    phi_idx = np.clip(
        np.floor(phi / np.pi * N_ANGLES).astype(int), 0, N_ANGLES - 1)
    r_idx   = np.floor((r   - r_edges[0]) / dr             ).astype(int)
    z_idx   = np.floor((z_mid - z_edges[0]) / slice_thickness).astype(int)
 
    in_bounds = (valid
                 & (r_idx >= 0) & (r_idx < N_RADIAL)
                 & (z_idx >= 0) & (z_idx < Z_RES))
 
    sinograms = np.zeros((Z_RES, N_ANGLES, N_RADIAL), dtype=np.float64)
    np.add.at(sinograms,
              (z_idx[in_bounds], phi_idx[in_bounds], r_idx[in_bounds]),
              1.0)
 
    print(f"Binned {in_bounds.sum():,} LORs into "
          f"{Z_RES} x {N_ANGLES} x {N_RADIAL} sinogram")
    return sinograms
 
 
# ============================================================
# Step 4: Ramp filter (Ram-Lak) with Hann apodization
# ============================================================
def ramp_filter(sinogram_slice):
    """
    Ram-Lak ramp filter with Hann apodization along radial axis.
    Hann window suppresses high-frequency streak artifacts from
    non-uniform angular sampling without significantly broadening the PSF.
    """
    S     = np.fft.rfft(sinogram_slice, axis=1)
    freqs = np.fft.rfftfreq(N_RADIAL, d=dr)
    ramp  = np.abs(freqs)
    n     = len(freqs)
    hann  = 0.5 * (1.0 - np.cos(np.pi * np.arange(n) / (n - 1)))
    return np.fft.irfft(S * ramp * hann[np.newaxis, :],
                        n=N_RADIAL, axis=1)
 
 
# ============================================================
# Step 5: Backproject one filtered sinogram into a 2D image
# ============================================================
def backproject(filtered_sinogram):
    """
    2D filtered backprojection for one axial slice.
    Uses the same angle convention as build_sinograms:
        phi = arctan2(lx, -ly)  →  cos(phi) = -ly, sin(phi) = lx
        r   = x*cos(phi) + y*sin(phi) = x*(-ly) + y*lx = x1*ly - y1*lx  ✓
    """
    image = np.zeros((X_RES, Y_RES), dtype=np.float64)
    for a_idx in range(N_ANGLES):
        r_pixel = PX * COS_PHI[a_idx] + PY * SIN_PHI[a_idx]
        r_frac  = (r_pixel - r_edges[0]) / dr
        r_lo    = np.floor(r_frac).astype(int)
        r_hi    = r_lo + 1
        w_hi    = r_frac - r_lo
        w_lo    = 1.0 - w_hi
        proj    = filtered_sinogram[a_idx]
        in_b    = (r_lo >= 0) & (r_hi < N_RADIAL)
        image  += np.where(
            in_b,
            w_lo * proj[np.clip(r_lo, 0, N_RADIAL - 1)]
          + w_hi * proj[np.clip(r_hi, 0, N_RADIAL - 1)],
            0.0)
    image *= np.pi / N_ANGLES
    return image
 
 
# ============================================================
# Step 6: Write .voxels file
# ============================================================
def write_voxels(filepath, volume):
    with open(filepath, 'wb') as f:
        np.array([X_LEN], dtype=np.float64).tofile(f)
        np.array([Y_LEN], dtype=np.float64).tofile(f)
        np.array([Z_LEN], dtype=np.float64).tofile(f)
        np.array([1],     dtype=np.int32  ).tofile(f)  # 1 iteration for FBP
        np.array([X_RES], dtype=np.int32  ).tofile(f)
        np.array([Y_RES], dtype=np.int32  ).tofile(f)
        np.array([Z_RES], dtype=np.int32  ).tofile(f)
        volume.astype(np.float64).tofile(f)
    print(f"Written: {filepath}")
    print(f"  Volume: {X_RES} x {Y_RES} x {Z_RES} voxels")
    print(f"  FOV:    {X_LEN:.1f} x {Y_LEN:.1f} x {Z_LEN:.1f} cm")
    print(f"  Voxel:  {VOXEL_SIZE_CM * 10:.1f} mm isotropic")
 
 
# ============================================================
# Main
# ============================================================
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python fbp_reconstruction.py "
              "<lor_endpoint_file> <output.voxels>")
        sys.exit(1)
 
    lor_file    = sys.argv[1]
    output_file = sys.argv[2]
 
    print("=" * 55)
    print("FBP Reconstruction")
    print("=" * 55)
    print(f"Reconstruction volume: "
          f"{X_LEN:.1f} x {Y_LEN:.1f} x {Z_LEN:.1f} cm")
    print(f"Voxel size:  {VOXEL_SIZE_CM * 10:.1f} mm")
    print(f"Matrix size: {X_RES} x {Y_RES} x {Z_RES}")
    print(f"Sinogram:    {N_ANGLES} angles x {N_RADIAL} radial bins")
    print(f"Axial cut:   {MAX_AXIAL_ANGLE_DEG:.1f} degrees")
    print()
 
    lors = read_lors(lor_file)
    lors = axial_cut(lors)
 
    if len(lors) == 0:
        print("ERROR: No LORs passed the axial cut.")
        print("Try increasing MAX_AXIAL_ANGLE_DEG and retry.")
        sys.exit(1)
 
    print("\nBuilding sinograms...")
    sinograms = build_sinograms(lors)
 
    print("Filtering and backprojecting...")
    volume = np.zeros((X_RES, Y_RES, Z_RES), dtype=np.float64)
    for z in range(Z_RES):
        if z % 20 == 0:
            print(f"  Slice {z + 1}/{Z_RES}")
        if sinograms[z].sum() == 0:
            continue
        filtered        = ramp_filter(sinograms[z])
        volume[:, :, z] = backproject(filtered)
 
    np.clip(volume, 0, None, out=volume)
 
    print("\nWriting output...")
    write_voxels(output_file, volume)
    print("\nDone. Visualize with:")
    print(f"  python render.py {output_file}")
 