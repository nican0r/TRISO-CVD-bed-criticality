"""Visualise how the tiling bounding-box grid covers the cone surface.

Geometry (from params.yaml, converted to cm):
  r_throat  = 0.3 cm        (z = 0)
  r_retort  = 2.5 cm        (z = z_cone_top)
  z_cone_top = 3.8 cm
  half-angle = 30°  →  tan θ = 1/√3

Tiling lattice (tile_size_cm = 0.5, used in exact_cone_bed):
  lower_left  = (−r_retort, −r_retort, 0)
  upper_right = ( r_retort,  r_retort, z_bed_top)
  nx = ny = ⌈5.0 / 0.5⌉ = 10,   nz = ⌈3.8 / 0.5⌉ = 8

Tile status per cell:
  fully inside   — all 8 corners inside the frustum
  boundary (wall)— some corners inside, some outside
  fully outside  — all corners outside the frustum
"""

import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d.proj3d import proj_transform
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

matplotlib.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 9,
    'axes.titlesize': 10,
    'axes.labelsize': 9,
})

# ---------------------------------------------------------------------------
# Geometry constants (params.yaml → cm)
# ---------------------------------------------------------------------------
R_THROAT   = 0.3       # cm  (r at z=0)
R_RETORT   = 2.5       # cm  (r at z_cone_top)
Z_CONE_TOP = 3.8       # cm
HALF_ANGLE = 30.0      # degrees
TAN_THETA  = math.tan(math.radians(HALF_ANGLE))
Z_APEX     = -R_THROAT / TAN_THETA          # cm  (below z=0)

TILE_SIZE  = 0.5       # cm

# Bounding box of the tiling lattice
BB_X_LO, BB_X_HI = -R_RETORT, R_RETORT     # –2.5 … 2.5 cm
BB_Y_LO, BB_Y_HI = -R_RETORT, R_RETORT
BB_Z_LO           = 0.0
NX = NY = math.ceil((BB_X_HI - BB_X_LO) / TILE_SIZE)   # 10
NZ = math.ceil((Z_CONE_TOP - BB_Z_LO)   / TILE_SIZE)    # 8
BB_Z_HI = BB_Z_LO + NZ * TILE_SIZE                       # 4.0 cm


def cone_r(z):
    """Inner cone radius at height z."""
    return TAN_THETA * (z - Z_APEX)


def point_inside_frustum(x, y, z):
    """True if (x,y,z) is inside the truncated cone (frustum) region.

    Region: –cone_surface & +z_bot & –z_top  (OpenMC sign convention)
    """
    if z < 0.0 or z > Z_CONE_TOP:
        return False
    rho = math.hypot(x, y)
    return rho <= cone_r(z)


# ---------------------------------------------------------------------------
# Classify every tile
# ---------------------------------------------------------------------------
# Corners of a tile at index (ix, iy, iz):
#   x ∈ [BB_X_LO + ix*T,   BB_X_LO + (ix+1)*T]
#   z ∈ [BB_Z_LO + iz*T,   BB_Z_LO + (iz+1)*T]

def tile_status(ix, iy, iz):
    """Return 'inside', 'boundary', or 'outside' for a 3-D tile."""
    T = TILE_SIZE
    xs = [BB_X_LO + ix * T, BB_X_LO + (ix + 1) * T]
    ys = [BB_Y_LO + iy * T, BB_Y_LO + (iy + 1) * T]
    zs = [BB_Z_LO + iz * T, BB_Z_LO + (iz + 1) * T]
    corners = [
        point_inside_frustum(x, y, z)
        for x in xs for y in ys for z in zs
    ]
    n_in = sum(corners)
    if n_in == 8:
        return 'inside'
    elif n_in == 0:
        return 'outside'
    else:
        return 'boundary'


# Build status grid
status = np.empty((NX, NY, NZ), dtype=object)
for ix in range(NX):
    for iy in range(NY):
        for iz in range(NZ):
            status[ix, iy, iz] = tile_status(ix, iy, iz)

# ---------------------------------------------------------------------------
# Colour map
# ---------------------------------------------------------------------------
COLORS = {
    'inside':   '#4393c3',   # blue
    'boundary': '#f4a582',   # salmon / orange
    'outside':  '#f7f7f7',   # near-white
}
ALPHAS = {'inside': 0.55, 'boundary': 0.80, 'outside': 0.10}
EC     = {'inside': '#2166ac', 'boundary': '#d6604d', 'outside': '#cccccc'}

# ---------------------------------------------------------------------------
# Figure layout: two panels (2D cross-section + 3D perspective)
# ---------------------------------------------------------------------------
fig = plt.figure(figsize=(13, 6.5))
ax2d = fig.add_subplot(1, 2, 1)
ax3d = fig.add_subplot(1, 2, 2, projection='3d')

# ============================================================
# Panel A — 2D cross-section (x–z plane, y = 0 row of tiles)
# ============================================================
T = TILE_SIZE
z_arr = np.linspace(0.0, Z_CONE_TOP, 300)
r_arr = [cone_r(z) for z in z_arr]

iy_mid = NY // 2   # middle y row  (y ∈ [0.0, 0.5] cm when NY=10)

for ix in range(NX):
    for iz in range(NZ):
        st = status[ix, iy_mid, iz]
        x0 = BB_X_LO + ix * T
        z0 = BB_Z_LO + iz * T
        rect = plt.Rectangle(
            (x0, z0), T, T,
            facecolor=COLORS[st], edgecolor=EC[st],
            linewidth=0.6, alpha=ALPHAS[st],
        )
        ax2d.add_patch(rect)

# Cone walls
ax2d.plot( r_arr, z_arr, color='#333333', lw=2.0, label='Cone wall')
ax2d.plot([-v for v in r_arr], z_arr, color='#333333', lw=2.0)

# Bounding box outline
ax2d.add_patch(plt.Rectangle(
    (BB_X_LO, BB_Z_LO), BB_X_HI - BB_X_LO, BB_Z_HI - BB_Z_LO,
    fill=False, edgecolor='#555555', lw=1.2, linestyle='--', label='Tiling bounding box',
))

# Annotations for r_throat and r_retort
ax2d.annotate(
    '', xy=(R_THROAT, 0), xytext=(0, 0),
    arrowprops=dict(arrowstyle='<->', color='#666666', lw=1.0),
)
ax2d.text(R_THROAT / 2, -0.15, f'r_throat\n{R_THROAT} cm', ha='center',
          va='top', fontsize=7.5, color='#444444')

ax2d.annotate(
    '', xy=(R_RETORT, Z_CONE_TOP), xytext=(0, Z_CONE_TOP),
    arrowprops=dict(arrowstyle='<->', color='#666666', lw=1.0),
)
ax2d.text(R_RETORT / 2, Z_CONE_TOP + 0.18, f'r_retort = {R_RETORT} cm',
          ha='center', va='bottom', fontsize=7.5, color='#444444')

# Tile-grid lines (faint)
for i in range(NX + 1):
    xg = BB_X_LO + i * T
    ax2d.axvline(xg, color='#aaaaaa', lw=0.3, zorder=0)
for k in range(NZ + 1):
    zg = BB_Z_LO + k * T
    ax2d.axhline(zg, color='#aaaaaa', lw=0.3, zorder=0)

# Legend patches
patches = [
    mpatches.Patch(facecolor=COLORS['inside'],   edgecolor=EC['inside'],   alpha=0.75, label='Fully inside cone'),
    mpatches.Patch(facecolor=COLORS['boundary'], edgecolor=EC['boundary'], alpha=0.85, label='Straddles cone wall'),
    mpatches.Patch(facecolor=COLORS['outside'],  edgecolor=EC['outside'],  alpha=0.50, label='Outside cone (inert)'),
    mpatches.Patch(fill=False, edgecolor='#333333', lw=1.5, label='Cone wall'),
    mpatches.Patch(fill=False, edgecolor='#555555', linestyle='--', lw=1.2, label='Tiling bounding box'),
]
ax2d.legend(handles=patches, loc='upper left', fontsize=7.5, framealpha=0.9)

ax2d.set_xlim(BB_X_LO - 0.1, BB_X_HI + 0.1)
ax2d.set_ylim(-0.4, BB_Z_HI + 0.3)
ax2d.set_xlabel('x  (cm)')
ax2d.set_ylabel('z  (cm)')
ax2d.set_aspect('equal')
ax2d.set_title(
    f'Cross-section (x–z, y-centre row)\n'
    f'{NX}×{NZ} tiles visible, tile size = {T} cm',
)

n_in   = int((status == 'inside').sum())
n_bd   = int((status == 'boundary').sum())
n_out  = int((status == 'outside').sum())
n_tot  = NX * NY * NZ
ax2d.text(0.02, 0.03,
    f'Total tiles: {n_tot}  ({NX}×{NY}×{NZ})\n'
    f'Inside: {n_in}  |  Boundary: {n_bd}  |  Outside: {n_out}',
    transform=ax2d.transAxes, fontsize=7.5,
    va='bottom', ha='left',
    bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='#cccccc', alpha=0.9),
)

# ============================================================
# Panel B — 3D perspective
# ============================================================
ax3d.set_box_aspect([1, 1, 1.0])

def draw_tile_3d(ax, ix, iy, iz, fc, ec, alpha, lw=0.4):
    """Draw one cubic tile as six Poly3DCollection faces."""
    x0 = BB_X_LO + ix * T;  x1 = x0 + T
    y0 = BB_Y_LO + iy * T;  y1 = y0 + T
    z0 = BB_Z_LO + iz * T;  z1 = z0 + T
    verts = [
        [(x0,y0,z0),(x1,y0,z0),(x1,y1,z0),(x0,y1,z0)],  # bottom
        [(x0,y0,z1),(x1,y0,z1),(x1,y1,z1),(x0,y1,z1)],  # top
        [(x0,y0,z0),(x1,y0,z0),(x1,y0,z1),(x0,y0,z1)],  # front
        [(x0,y1,z0),(x1,y1,z0),(x1,y1,z1),(x0,y1,z1)],  # back
        [(x0,y0,z0),(x0,y1,z0),(x0,y1,z1),(x0,y0,z1)],  # left
        [(x1,y0,z0),(x1,y1,z0),(x1,y1,z1),(x1,y0,z1)],  # right
    ]
    poly = Poly3DCollection(verts, facecolor=fc, edgecolor=ec,
                            linewidth=lw, alpha=alpha)
    ax.add_collection3d(poly)


for ix in range(NX):
    for iy in range(NY):
        for iz in range(NZ):
            st = status[ix, iy, iz]
            if st == 'outside':
                continue   # don't draw outside tiles in 3D (too cluttered)
            draw_tile_3d(ax3d, ix, iy, iz,
                         fc=COLORS[st], ec=EC[st],
                         alpha=ALPHAS[st] + 0.1,
                         lw=0.3)

# Draw cone surface as a wire-frame frustum
theta_vals = np.linspace(0, 2 * math.pi, 60)
z_vals     = np.linspace(0, Z_CONE_TOP, 20)
for z in z_vals:
    r = cone_r(z)
    xs = r * np.cos(theta_vals)
    ys = r * np.sin(theta_vals)
    ax3d.plot(xs, ys, z, color='#222222', lw=0.7, alpha=0.6)
for theta in np.linspace(0, 2 * math.pi, 12, endpoint=False):
    rs = [cone_r(z) for z in z_vals]
    xs = [r * math.cos(theta) for r in rs]
    ys = [r * math.sin(theta) for r in rs]
    ax3d.plot(xs, ys, z_vals, color='#222222', lw=0.7, alpha=0.6)

# Bounding-box outline (dashed)
corners_bb = [
    (BB_X_LO, BB_Y_LO, BB_Z_LO), (BB_X_HI, BB_Y_LO, BB_Z_LO),
    (BB_X_HI, BB_Y_HI, BB_Z_LO), (BB_X_LO, BB_Y_HI, BB_Z_LO),
    (BB_X_LO, BB_Y_LO, BB_Z_LO),
]
for a, b in zip(corners_bb, corners_bb[1:]):
    ax3d.plot([a[0],b[0]], [a[1],b[1]], [a[2],b[2]], '--', color='#888888', lw=0.8)
corners_top = [
    (BB_X_LO, BB_Y_LO, BB_Z_HI), (BB_X_HI, BB_Y_LO, BB_Z_HI),
    (BB_X_HI, BB_Y_HI, BB_Z_HI), (BB_X_LO, BB_Y_HI, BB_Z_HI),
    (BB_X_LO, BB_Y_LO, BB_Z_HI),
]
for a, b in zip(corners_top, corners_top[1:]):
    ax3d.plot([a[0],b[0]], [a[1],b[1]], [a[2],b[2]], '--', color='#888888', lw=0.8)
for c_b, c_t in zip(corners_bb[:4], corners_top[:4]):
    ax3d.plot([c_b[0],c_t[0]], [c_b[1],c_t[1]], [c_b[2],c_t[2]], '--',
              color='#888888', lw=0.8)

ax3d.set_xlabel('x  (cm)', labelpad=4)
ax3d.set_ylabel('y  (cm)', labelpad=4)
ax3d.set_zlabel('z  (cm)', labelpad=4)
ax3d.set_title(
    f'3-D tiling grid (inside + boundary tiles only)\n'
    f'{n_in} inside  +  {n_bd} boundary tiles drawn',
)
ax3d.view_init(elev=25, azim=-50)
ax3d.set_xlim(BB_X_LO, BB_X_HI)
ax3d.set_ylim(BB_Y_LO, BB_Y_HI)
ax3d.set_zlim(BB_Z_LO, BB_Z_HI)

fig.suptitle(
    'Tiling bounding-box coverage of the CVD furnace cone\n'
    f'Tile size: {T} cm × {T} cm × {T} cm   |   '
    f'Grid: {NX}×{NY}×{NZ} = {n_tot} tiles   |   '
    f'Cone: 60° full angle, {R_THROAT}–{R_RETORT} cm radius, {Z_CONE_TOP} cm tall',
    fontsize=9.5, y=1.01,
)

plt.tight_layout()
out_path = 'plots/tiling_cone_coverage.png'
import pathlib; pathlib.Path('plots').mkdir(exist_ok=True)
plt.savefig(out_path, dpi=160, bbox_inches='tight')
print(f"Saved → {out_path}")
plt.show()
