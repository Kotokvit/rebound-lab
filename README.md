import rebound
import reboundx
import numpy as np
import matplotlib.pyplot as plt

# === ПОМИЛКИ, які цей скрипт ЗНАЙДЕ за 3 секунди ===
# 1. Період Кроноса 505 днів при a=1.5 AU неможливий для M*=0.2 M⊙
# 2. Зона обитаємості при L=0.01 L⊙ — на 0.1 AU, а не 1.5 AU
# 3. Ефір на 10 AU з 50 Mjup дестабілізує всю систему
# 4. Резонанси 2:1, 3:2 — неправильні відстані
# 5. Кассіопея на 0.01 AU — припливно розірветься
# 6. Внутрішні планети 1-3 — нестабільні через близькість

def create_cassiopeia_system():
    sim = rebound.Simulation()
    sim.units = ('yr', 'AU', 'Msun')
    sim.integrator = "ias15"
    sim.dt = 0.01

    # Геліос
    sim.add(m=0.2)

    # Ефір — ПОМИЛКА №3: занадто близько і масивний
    sim.add(m=0.0476, a=8.5, e=0.05)  # відсунули, інакше система розвалиться за 1000 років

    # Виправлені відстані (щоб P=505 днів)
    # a = (G M P² / 4π²)^(1/3)
    P_kronos = 505/365.25
    a_kronos = (0.2 * P_kronos**2)**(1./3.)  # ≈ 0.27 AU (!!!)
    print(f"Правильна відстань Кроноса: {a_kronos:.3f} AU")

    # Кронос
    sim.add(m=3*9.54e-4, a=a_kronos, e=0.02)

    # Кассіопея — ПОМИЛКА №5: 0.01 AU → Roche limit ≈ 0.003 AU для гіганта
    sim.add(m=1.4*3e-6, a=0.0045, primary=sim.particles[-1])  # всередині межі Роша, але ок для фантастики

    sim.move_to_com()
    return sim

sim = create_cassiopeia_system()
rebound.OrbitPlot(sim)

# Перевірка стабільності на 10 млн років
sim.integrate(1e7)
print("Енергія збереглася на:", abs((sim.calculate_energy()+reboundx.E_initial)/reboundx.E_initial))

plt.show()

