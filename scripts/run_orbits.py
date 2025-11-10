#!/usr/bin/env python3
import argparse, os, math, json
import numpy as np, pandas as pd, matplotlib.pyplot as plt
import rebound

def add_system(sim):
    Mstar = 0.20  # Msun
    sim.add(m=Mstar)
    def MJ(x): return x*0.0009543
    sim.add(m=MJ(0.05), a=0.63, e=0.02, hash="p6")
    sim.add(m=MJ(3.0),  a=0.73, e=0.03, hash="kronos")
    sim.add(m=MJ(0.15), a=0.90, e=0.02, hash="p8")
    sim.add(m=MJ(0.12), a=1.43, e=0.02, hash="p9")
    return Mstar

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--years", type=float, default=5000)
    ap.add_argument("--dt", type=float, default=5e-4)
    ap.add_argument("--out", type=str, default="artifacts")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    sim = rebound.Simulation()
    sim.units = ('AU','yr','Msun')
    Mstar = add_system(sim)
    sim.move_to_com()
    sim.integrator = "whfast"
    sim.dt = args.dt

    ts = np.linspace(0, args.years, int(args.years/args.dt/50)+1)
    rows = []
    for t in ts:
        sim.integrate(t)
        o6 = rebound.Orbit(sim.particles["p6"], sim.particles[0])
        o7 = rebound.Orbit(sim.particles["kronos"], sim.particles[0])
        o8 = rebound.Orbit(sim.particles["p8"], sim.particles[0])
        o9 = rebound.Orbit(sim.particles["p9"], sim.particles[0])
        phi_78 = (3*o8.l - 2*o7.l - o8.pomega)%(2*np.pi)
        phi_67 = (2*o7.l - 1*o6.l - o7.pomega)%(2*np.pi)
        rows.append(dict(t=t,a6=o6.a,e6=o6.e,a7=o7.a,e7=o7.e,a8=o8.a,e8=o8.e,a9=o9.a,e9=o9.e,
                         phi_78=phi_78,phi_67=phi_67,energy=sim.calculate_energy()))
    import pandas as pd
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(args.out,"orbits.csv"), index=False)
    plt.figure(); 
    for lab in ["a6","a7","a8","a9"]: plt.plot(df["t"], df[lab], label=lab)
    plt.xlabel("Time [yr]"); plt.ylabel("a [AU]"); plt.legend(); plt.tight_layout()
    plt.savefig(os.path.join(args.out,"a_vs_time.png"), dpi=200)
    plt.figure(); 
    for lab in ["e6","e7","e8","e9"]: plt.plot(df["t"], df[lab], label=lab)
    plt.xlabel("Time [yr]"); plt.ylabel("e"); plt.legend(); plt.tight_layout()
    plt.savefig(os.path.join(args.out,"e_vs_time.png"), dpi=200)
    E0 = df["energy"].iloc[0]; dE = (df["energy"]-E0)/abs(E0)
    plt.figure(); plt.plot(df["t"], dE); plt.xlabel("Time [yr]"); plt.ylabel("Rel energy drift"); plt.tight_layout()
    plt.savefig(os.path.join(args.out,"energy_drift.png"), dpi=200)
    with open(os.path.join(args.out,"summary.json"),"w",encoding="utf-8") as f:
        import json; json.dump(dict(energy_rel_drift=float(abs(dE).max())), f, indent=2)
    print("[OK] Wrote artifacts to", args.out)

if __name__ == "__main__":
    main()
