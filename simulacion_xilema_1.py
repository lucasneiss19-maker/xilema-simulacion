"""
Simulación: Capilaridad y Flujo en Xilema
==========================================
Autores: Lucas Neisskenwirth, Martín Costa, Max Lavandero
Universidad Adolfo Ibáñez, 2026

Modelo de red de poros con Hagen-Poiseuille.
Genera todos los gráficos y la tabla de resultados en /outputs/.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.linalg import solve
import warnings
warnings.filterwarnings("ignore")

np.random.seed(42)  # reproducibilidad

# ──────────────────────────────────────────────
# PARÁMETROS FÍSICOS (SI)
# ──────────────────────────────────────────────
GAMMA   = 0.0728   # N/m  — tensión superficial agua a 20°C
THETA   = 0.0      # rad  — ángulo de contacto (cos θ ≈ 1)
RHO     = 998.0    # kg/m³ — densidad agua
G       = 9.81     # m/s²
ETA     = 1e-3     # Pa·s — viscosidad dinámica agua a 20°C
L_VASO  = 0.01     # m    — longitud de cada segmento de vaso (1 cm)

# Radios de vasos (m) — rango típico xilema
R_TRACHEID   = 10e-6   # 10 µm
R_VESSEL_TYP = 20e-6   # 20 µm (radio de referencia para Jurin)
R_VESSEL_LRG = 50e-6   # 50 µm
R_MAX        = 500e-6  # límite superior para barrido

# ──────────────────────────────────────────────
# 1. LEY DE JURIN: altura capilar vs radio
# ──────────────────────────────────────────────

def jurin_height(r):
    """Altura capilar (m) para radio r (m)."""
    return (2 * GAMMA * np.cos(THETA)) / (RHO * G * r)


radii_jurin = np.linspace(1e-6, 200e-6, 500)  # 1–200 µm
heights_m   = jurin_height(radii_jurin)

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(radii_jurin * 1e6, heights_m * 100, color="#2E7D52", lw=2.5)
ax.axvline(20, color="#B03A2E", ls="--", lw=1.5, label="Radio típico xilema (20 µm)")
ax.axhline(jurin_height(20e-6) * 100, color="#B03A2E", ls=":", lw=1.2)
ax.annotate(f"h ≈ {jurin_height(20e-6)*100:.1f} cm\n(r = 20 µm)",
            xy=(20, jurin_height(20e-6)*100),
            xytext=(60, jurin_height(20e-6)*100 + 30),
            fontsize=9,
            arrowprops=dict(arrowstyle="->", color="#B03A2E"),
            color="#B03A2E")
ax.set_xlabel("Radio del vaso (µm)", fontsize=12)
ax.set_ylabel("Altura capilar (cm)", fontsize=12)
ax.set_title("Ley de Jurin: Ascenso capilar en xilema", fontsize=13, fontweight="bold")
ax.legend(fontsize=9)
ax.set_xlim(0, 200)
ax.set_ylim(0, heights_m.max() * 100 * 1.05)
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("../outputs/fig1_jurin_ascenso.png", dpi=150)
plt.close()
print("✓ fig1_jurin_ascenso.png")

# ──────────────────────────────────────────────
# 2. HAGEN-POISEUILLE: caudal vs radio
# ──────────────────────────────────────────────

def hagen_poiseuille(r, dP, L=L_VASO):
    """Caudal Q (m³/s) según Hagen-Poiseuille."""
    return (np.pi * r**4 * dP) / (8 * ETA * L)


DELTA_P = 1e6   # Pa (≈ 10 bar — gradiente típico xilema)
radii_hp = np.linspace(1e-6, 100e-6, 500)
Q_vals   = hagen_poiseuille(radii_hp, DELTA_P)

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(radii_hp * 1e6, Q_vals * 1e12, color="#1A5276", lw=2.5)
# Marcar r doble → Q×16
r_ref = 25e-6
Q_ref = hagen_poiseuille(r_ref, DELTA_P)
Q_2r  = hagen_poiseuille(2*r_ref, DELTA_P)
ax.annotate("", xy=(2*r_ref*1e6, Q_2r*1e12),
            xytext=(r_ref*1e6, Q_ref*1e12),
            arrowprops=dict(arrowstyle="->", color="#884EA0", lw=2))
ax.text(2*r_ref*1e6 + 2, Q_2r*1e12, f"×16 en Q\nal doblar r", fontsize=9, color="#884EA0")
ax.set_xlabel("Radio del vaso (µm)", fontsize=12)
ax.set_ylabel("Caudal Q (pL/s)", fontsize=12)
ax.set_title("Hagen-Poiseuille: Caudal vs radio en xilema\n"
             r"Q = $\pi r^4 \Delta P / (8 \eta L)$  — $\Delta P$ = 1 MPa", fontsize=11, fontweight="bold")
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("../outputs/fig2_hagen_poiseuille.png", dpi=150)
plt.close()
print("✓ fig2_hagen_poiseuille.png")

# ──────────────────────────────────────────────
# 3. RED DE POROS: modelo compartimental
#    Nodos en columna vertical con condiciones
#    de contorno en raíz (nodo 0) y hoja (nodo N)
# ──────────────────────────────────────────────

def simular_red(N_nodos=10, radii_list=None, psi_raiz=0.0, psi_hoja=-1e6):
    """
    Red lineal de N_nodos segmentos de xilema.
    Resuelve Kirchhoff para potencial hídrico en cada nodo.
    Retorna nodos, potenciales y caudales.
    """
    if radii_list is None:
        radii_list = [R_VESSEL_TYP] * N_nodos

    N = len(radii_list)
    # Conductancias: K_i = pi*r^4 / (8*eta*L)
    K = np.array([np.pi * r**4 / (8 * ETA * L_VASO) for r in radii_list])

    # Ensamble matriz de Kirchhoff (nodos internos)
    n_int = N - 1  # nodos internos (0 y N son contorno)
    A = np.zeros((n_int, n_int))
    b = np.zeros(n_int)

    for i in range(n_int):
        A[i, i] += K[i] + K[i+1]
        if i > 0:
            A[i, i-1] -= K[i]
        if i < n_int - 1:
            A[i, i+1] -= K[i+1]
        if i == 0:
            b[i] += K[i] * psi_raiz
        if i == n_int - 1:
            b[i] += K[i+1] * psi_hoja

    psi_int = solve(A, b)
    psi_all = np.concatenate([[psi_raiz], psi_int, [psi_hoja]])

    # Caudal en cada segmento
    Q_segs = K * np.diff(psi_all)
    return np.arange(N+1), psi_all, Q_segs


# ── 3a. Escenario base: vasos uniformes 20 µm ──
nodos, psi, Q_segs = simular_red(N_nodos=20, psi_raiz=0.0, psi_hoja=-1e6)
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

ax = axes[0]
ax.plot(nodos * L_VASO * 100, psi / 1e6, "o-", color="#2E7D52", lw=2, ms=5)
ax.set_xlabel("Posición en tallo (cm)", fontsize=12)
ax.set_ylabel("Potencial hídrico ψ (MPa)", fontsize=12)
ax.set_title("Gradiente de potencial hídrico\n(red uniforme, r = 20 µm)", fontsize=11, fontweight="bold")
ax.grid(alpha=0.3)

ax = axes[1]
pos_mid = (nodos[:-1] + nodos[1:]) / 2 * L_VASO * 100
ax.bar(pos_mid, np.abs(Q_segs) * 1e12, width=0.7, color="#1A5276", alpha=0.8)
ax.set_xlabel("Posición en tallo (cm)", fontsize=12)
ax.set_ylabel("Caudal |Q| (pL/s)", fontsize=12)
ax.set_title("Caudal por segmento\n(estado estacionario)", fontsize=11, fontweight="bold")
ax.grid(alpha=0.3, axis="y")

plt.suptitle("Red de poros — Modelo base (20 segmentos, ΔΨ = 1 MPa)", fontsize=12)
plt.tight_layout()
plt.savefig("../outputs/fig3_red_poros_base.png", dpi=150)
plt.close()
print("✓ fig3_red_poros_base.png")

# ──────────────────────────────────────────────
# 4. ANÁLISIS DE SENSIBILIDAD: barrido r y ΔΨ
# ──────────────────────────────────────────────

radii_sweep = np.array([5, 10, 20, 50, 100, 200, 500]) * 1e-6  # µm
psi_sweep   = np.array([0.1, 0.5, 1.0, 2.0, 5.0, 10.0]) * 1e6  # MPa

Q_matrix = np.zeros((len(psi_sweep), len(radii_sweep)))
for i, dp in enumerate(psi_sweep):
    for j, r in enumerate(radii_sweep):
        _, _, Q_s = simular_red(N_nodos=10,
                                radii_list=[r]*10,
                                psi_raiz=0.0,
                                psi_hoja=-dp)
        Q_matrix[i, j] = np.abs(Q_s).mean() * 1e12  # pL/s

fig, ax = plt.subplots(figsize=(9, 6))
im = ax.imshow(Q_matrix, aspect="auto", cmap="YlOrRd", origin="lower")
ax.set_xticks(range(len(radii_sweep)))
ax.set_xticklabels([f"{int(r*1e6)}" for r in radii_sweep])
ax.set_yticks(range(len(psi_sweep)))
ax.set_yticklabels([f"{p/1e6:.1f}" for p in psi_sweep])
ax.set_xlabel("Radio del vaso (µm)", fontsize=12)
ax.set_ylabel("Diferencia de potencial hídrico ΔΨ (MPa)", fontsize=12)
ax.set_title("Mapa de calor: Caudal medio (pL/s)\nen función de r y ΔΨ", fontsize=12, fontweight="bold")
for i in range(len(psi_sweep)):
    for j in range(len(radii_sweep)):
        val = Q_matrix[i, j]
        txt = f"{val:.1f}" if val < 1000 else f"{val/1000:.1f}k"
        ax.text(j, i, txt, ha="center", va="center",
                fontsize=7, color="black" if val < Q_matrix.max()*0.6 else "white")
plt.colorbar(im, ax=ax, label="Caudal medio (pL/s)")
plt.tight_layout()
plt.savefig("../outputs/fig4_sensibilidad_heatmap.png", dpi=150)
plt.close()
print("✓ fig4_sensibilidad_heatmap.png")

# ──────────────────────────────────────────────
# 5. COMPARACIÓN: traqueidas vs elementos de vaso
# ──────────────────────────────────────────────

configs = {
    "Solo traqueidas (10 µm)":      [R_TRACHEID]   * 20,
    "Vasos típicos (20 µm)":        [R_VESSEL_TYP] * 20,
    "Vasos grandes (50 µm)":        [R_VESSEL_LRG] * 20,
    "Mixto (traqueidas+vasos)":     [R_TRACHEID if i % 2 == 0 else R_VESSEL_TYP for i in range(20)],
}

colors  = ["#884EA0", "#2E7D52", "#1A5276", "#B7950B"]
fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True)
axes = axes.flatten()

for ax, (label, rlist), col in zip(axes, configs.items(), colors):
    nds, psi_c, Q_c = simular_red(N_nodos=20, radii_list=rlist, psi_raiz=0.0, psi_hoja=-1e6)
    ax.plot(nds * L_VASO * 100, psi_c / 1e6, "o-", color=col, lw=2, ms=4)
    ax.set_title(label, fontsize=10, fontweight="bold")
    ax.set_ylabel("ψ (MPa)", fontsize=9)
    ax.set_xlabel("Posición (cm)", fontsize=9)
    ax.grid(alpha=0.3)
    ax.text(0.02, 0.05, f"Q̄ = {np.abs(Q_c).mean()*1e12:.1f} pL/s",
            transform=ax.transAxes, fontsize=8, color=col,
            bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.7))

plt.suptitle("Comparación de configuraciones anatómicas del xilema\n(ΔΨ = 1 MPa, 20 segmentos)",
             fontsize=12, fontweight="bold")
plt.tight_layout()
plt.savefig("../outputs/fig5_comparacion_configs.png", dpi=150)
plt.close()
print("✓ fig5_comparacion_configs.png")

# ──────────────────────────────────────────────
# 6. TABLA RESUMEN (CSV)
# ──────────────────────────────────────────────

import csv

rows = []
for label, rlist in configs.items():
    _, _, Q_c = simular_red(N_nodos=20, radii_list=rlist, psi_raiz=0.0, psi_hoja=-1e6)
    r_mean = np.mean(rlist) * 1e6
    Q_mean = np.abs(Q_c).mean() * 1e12
    h_cap  = jurin_height(np.mean(rlist)) * 100
    rows.append([label, f"{r_mean:.0f}", f"{Q_mean:.2f}", f"{h_cap:.1f}"])

with open("../outputs/tabla_resultados.csv", "w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f)
    writer.writerow(["Configuración", "Radio medio (µm)", "Caudal medio (pL/s)", "Altura Jurin (cm)"])
    writer.writerows(rows)

print("✓ tabla_resultados.csv")
print("\n✅ Simulación completa. Archivos en /outputs/")
