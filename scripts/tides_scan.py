#!/usr/bin/env python3
import argparse, os, numpy as np, pandas as pd, matplotlib.pyplot as plt

AU = 1.495978707e11
R_E = 6.371e6
M_J = 1.89813e27

def hill_radius(a_AU, Mp_MJ, Mstar_Msun):
    return a_AU * ((Mp_MJ*0.0009543)/(3*Mstar_Msun))**(1/3)

def tidal_flux(Mp_MJ=3.0, Rs_Re=1.2, a_s_AU=0.01, e=0.03, k2_over_Q=2e-3):
    G = 6.67430e-11; Rs = Rs_Re*R_E; Mp = Mp_MJ*M_J; a = a_s_AU*AU
    n = (G*Mp/a**3)**0.5
    Edot = 21/2 * k2_over_Q * G*Mp**2 * Rs**5 * n * e**2 / a**6
    return Edot/(4*np.pi*Rs**2)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--a_AU", type=float, default=0.73)
    ap.add_argument("--Mp_MJ", type=float, default=3.0)
    ap.add_argument("--Mstar", type=float, default=0.20)
    ap.add_argument("--out", type=str, default="artifacts")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    RH = hill_radius(args.a_AU, args.Mp_MJ, args.Mstar)
    a_fracs = np.linspace(0.02, 0.30, 40)
    es = np.linspace(0.0, 0.10, 41)
    grid = np.zeros((len(es), len(a_fracs)))
    for i,e in enumerate(es):
        for j,f in enumerate(a_fracs):
            a_s = f*RH
            grid[i,j] = tidal_flux(args.Mp_MJ, 1.2, a_s, e, 2e-3)
    import pandas as pd
    df = pd.DataFrame(grid, index=np.round(es,3), columns=np.round(a_fracs,3))
    df.to_csv(os.path.join(args.out,"tides_grid.csv"))
    plt.figure(); plt.imshow(grid, origin="lower", aspect="auto",
                             extent=[a_fracs.min(), a_fracs.max(), es.min(), es.max()])
    plt.colorbar(label="Flux [W/m^2]"); plt.xlabel("a_s/R_H"); plt.ylabel("e")
    plt.tight_layout(); plt.savefig(os.path.join(args.out,"tides_heatmap.png"), dpi=200)
    print("[OK] Tidal grid written.")

if __name__ == "__main__":
    main()
