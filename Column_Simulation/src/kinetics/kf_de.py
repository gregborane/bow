import multiprocessing

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import differential_evolution
from scipy.sparse import diags

# =============================================================================
# 1. PARAMÈTRES FIXES DE LA COLONNE ET DU PROCESSUS
# =============================================================================
L = 0.042  # Longueur de la colonne (m)
R = 45e-6  # Rayon moyen de la bille (m)
eps_c = 0.33  # Porosité externe/interstitielle (-)
eps_p = 0.48  # Porosité apparente des billes pour l'IgG (-)
v_int = 7.4e-4  # Vitesse interstitielle (m/s)
D_ax = 2.0e-7  # Coefficient de dispersion axiale (m^2/s)

a = 3.0 / R  # Surface spécifique externe
v = v_int * eps_c  # Vitesse superficielle

Nz = 30  # Discrétisation z
Nr = 10  # Discrétisation r

dz = L / (Nz - 1)
dr = R / (Nr - 1)

z_nodes = np.linspace(0, L, Nz)
r_nodes = np.linspace(0, R, Nr)
r_inner = r_nodes[1:-1]
N_tot = Nz + Nz * Nr

sparsity_structure = diags(
    [-1.0, 0.0, 1.0], [-1, 0, 1], shape=(N_tot, N_tot), format="csc", dtype=float
)


# =============================================================================
# 2. SYSTÈME D'ÉQUATIONS (Vectorisé)
# =============================================================================
def column_pde_system(t, y, k_f, D_e):
    c = y[:Nz]
    # On lit les données, mais on NE DOIT JAMAIS modifier c_p en place
    c_p = y[Nz:].reshape((Nz, Nr))

    dc_dt = np.zeros(Nz)
    dcp_dt = np.zeros((Nz, Nr))

    # --- A. Phase stationnaire : Diffusion sphérique ---
    # Centre (r = 0)
    dcp_dt[:, 0] = (D_e / eps_p) * (6 * (c_p[:, 1] - c_p[:, 0]) / dr**2)

    # Cœur (0 < r < R)
    d2cp_dr2 = (c_p[:, 2:] - 2 * c_p[:, 1:-1] + c_p[:, :-2]) / dr**2
    dcp_dr = (c_p[:, 2:] - c_p[:, :-2]) / (2 * dr)
    dcp_dt[:, 1:-1] = (D_e / eps_p) * (d2cp_dr2 + (2.0 / r_inner) * dcp_dr)

    # Surface (r = R) - Résolue rigoureusement via un Point Fantôme
    c_surf = c_p[:, -1]
    # Point virtuel à l'extérieur de la bille pour respecter la condition de Robin
    c_ghost = c_p[:, -2] + (2 * dr * k_f / D_e) * (c - c_surf)

    d2cp_dr2_surf = (c_ghost - 2 * c_surf + c_p[:, -2]) / dr**2
    dcp_dr_surf = (k_f / D_e) * (c - c_surf)
    dcp_dt[:, -1] = (D_e / eps_p) * (d2cp_dr2_surf + (2.0 / R) * dcp_dr_surf)

    # --- B. Phase mobile : Convection, Dispersion et Film ---
    c_in = 1.0

    # z = 0 (Entrée) - Schéma Upwind pour la convection
    dc_dz_in = (c[0] - c_in) * (v / D_ax)
    d2c_dz2_in = (c[1] - 2 * c[0] + (c[0] - dc_dz_in * dz)) / dz**2
    dc_dz_0 = (c[0] - c_in) / dz  # Upwind stable

    dc_dt[0] = (
        D_ax * d2c_dz2_in
        - (v / eps_c) * dc_dz_0
        - ((1 - eps_c) / eps_c) * a * k_f * (c[0] - c_surf[0])
    )

    # 0 < z < L (Milieu)
    d2c_dz2 = (c[2:] - 2 * c[1:-1] + c[:-2]) / dz**2
    dc_dz = (c[1:-1] - c[:-2]) / dz  # SCHÉMA UPWIND

    dc_dt[1:-1] = (
        D_ax * d2c_dz2
        - (v / eps_c) * dc_dz
        - ((1 - eps_c) / eps_c) * a * k_f * (c[1:-1] - c_surf[1:-1])
    )

    # z = L (Sortie Neumann)
    d2c_dz2_out = (2 * c[-2] - 2 * c[-1]) / dz**2
    dc_dt[-1] = D_ax * d2c_dz2_out - ((1 - eps_c) / eps_c) * a * k_f * (
        c[-1] - c_surf[-1]
    )

    return np.concatenate([dc_dt, dcp_dt.flatten()])


def simulate_breakthrough(k_f, D_e, t_eval):
    y0 = np.zeros(N_tot)
    sol = solve_ivp(
        column_pde_system,
        t_span=(t_eval[0], t_eval[-1]),
        y0=y0,
        t_eval=t_eval,
        args=(k_f, D_e),
        method="BDF",
        # jac_sparsity=sparsity_structure,
        rtol=1e-5,
        atol=1e-8,
    )
    return sol.y[Nz - 1, :]


# =============================================================================
# 3. OPTIMISATION MULTI-COEURS (Somme des carrés des erreurs)
# =============================================================================
def sum_of_squares_error(params, t_eval, c_experimental):
    """L'optimiseur global minimise une valeur scalaire (SSE)."""
    k_f = 10 ** params[0]
    D_e = 10 ** params[1]

    try:
        c_sim = simulate_breakthrough(k_f, D_e, t_eval)
        # Pénalité sévère si le solveur échoue (ex: paramètres irréalistes)
        if len(c_sim) != len(c_experimental):
            return 1e6
        return np.sum((c_sim - c_experimental) ** 2)
    except:
        return 1e6


# =============================================================================
# 4. EXÉCUTION PRINCIPALE (Indispensable pour le multiprocessing)
# =============================================================================
if __name__ == "__main__":
    t_exp = np.linspace(0, 300, 50)
    # Génération des données factices
    c_exp = simulate_breakthrough(
        k_f=2.8e-6, D_e=4.1e-12, t_eval=t_exp
    ) + np.random.normal(0, 0.01, len(t_exp))

    # Limites de recherche pour l'espace des paramètres (log10)
    bounds = [(-8, -3), (-14, -9)]

    print(f"Démarrage de l'optimisation sur {multiprocessing.cpu_count()} cœurs...")

    # workers=-1 utilise tous les coeurs disponibles
    # updating='deferred' est requis lorsque l'on utilise un calcul parallèle
    result = differential_evolution(
        sum_of_squares_error,
        bounds,
        args=(t_exp, c_exp),
        workers=-1,  # <--- MULTIPROCESSING ICI
        updating="deferred",  # <--- REQUIS POUR WORKERS > 1
        disp=True,
        popsize=1,  # Taille de la population (réduire pour accélérer, augmenter pour la précision)
    )

    k_f_opt = 10 ** result.x[0]
    D_e_opt = 10 ** result.x[1]

    print("\n--- RÉSULTATS D'OPTIMISATION ---")
    print(f"k_f ajusté : {k_f_opt:.3e} m/s")
    print(f"D_e ajusté : {D_e_opt:.3e} m^2/s\n")

    # Affichage final
    t_sim_dense = np.linspace(0, 300, 200)
    c_sim_dense = simulate_breakthrough(k_f_opt, D_e_opt, t_sim_dense)

    print(c_sim_dense)

    # %%

    plt.figure(figsize=(8, 5))
    plt.plot(t_exp, c_exp, "kx", label="Données exp")
    plt.plot(
        t_sim_dense,
        c_sim_dense,
        "r-",
        label=f"Ajustement ($k_f$={k_f_opt:.2e}, $D_e$={D_e_opt:.2e})",
    )
    plt.xlabel("Temps [s]")
    plt.ylabel("C / Cin")
    plt.title("Courbe de percée - Optimisation Multi-cœurs")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig("tries.png")
# %%

