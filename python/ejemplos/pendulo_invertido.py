"""Péndulo invertido sostenido midiendo solo el ángulo.

Un regulador LQR calcula el par a partir del estado, y como la velocidad no
se mide la reconstruye un observador. El péndulo arranca desviado y a los 1.5 s
recibe un empujón.

    python ejemplos/pendulo_invertido.py                # abre la gráfica
    python ejemplos/pendulo_invertido.py pendulo.png    # la guarda en un archivo
"""

import math
import sys

import matplotlib

if len(sys.argv) > 1:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt

from setra import Observador, Pendulo, discretizar, ganancia_observador, lqr

DT = 2e-3
PAR_MAX = 0.5
MASA, LONGITUD, FRICCION, GRAVEDAD = 0.2, 0.5, 1e-3, 9.81
INICIO_EMPUJON, FIN_EMPUJON, PAR_EMPUJON = 1.5, 1.55, 0.4

pendulo = Pendulo(MASA, LONGITUD, friccion=FRICCION, gravedad=GRAVEDAD, angulo=math.pi + 0.2)

# Modelo lineal alrededor de la vertical hacia arriba; el estado es la
# desviación del ángulo y su velocidad.
inercia = MASA * LONGITUD**2
a, b = discretizar([[0.0, 1.0], [GRAVEDAD / LONGITUD, -FRICCION / inercia]], [[0.0], [1.0 / inercia]], DT)
c = [[1.0, 0.0]]
(k,) = lqr(a, b, [[10.0, 0.0], [0.0, 0.1]], [[1.0]])
observador = Observador(a, b, c, ganancia_observador(a, c, [[1.0, 0.0], [0.0, 100.0]], [[0.01]]))

tiempos, angulos, velocidades, estimadas, pares = [], [], [], [], []
for n in range(round(3.0 / DT)):
    t = n * DT
    desviacion, velocidad = observador.x
    par = max(-PAR_MAX, min(PAR_MAX, -(k[0] * desviacion + k[1] * velocidad)))
    observador.actualizar(par, pendulo.angulo - math.pi)

    tiempos.append(t)
    angulos.append(math.degrees(pendulo.angulo - math.pi))
    velocidades.append(pendulo.velocidad)
    estimadas.append(velocidad)
    pares.append(par)

    # El empujón llega al péndulo sin que el observador lo sepa.
    empujon = PAR_EMPUJON if INICIO_EMPUJON <= t < FIN_EMPUJON else 0.0
    pendulo.paso(par + empujon, DT, subpasos=4)

print(f"ganancia del regulador: {k[0]:.3f} N·m/rad, {k[1]:.3f} N·m·s/rad")
print(f"desviación máxima tras el empujón: {max(map(abs, angulos[round(INICIO_EMPUJON / DT):])):.2f}°")
print(f"desviación final: {abs(angulos[-1]):.4f}°")

figura, (arriba, medio, abajo) = plt.subplots(3, 1, sharex=True, figsize=(8, 7))
arriba.plot(tiempos, angulos)
arriba.set_ylabel("desviación (°)")
medio.plot(tiempos, velocidades, label="real")
medio.plot(tiempos, estimadas, "--", label="estimada")
medio.set_ylabel("velocidad (rad/s)")
medio.legend()
abajo.plot(tiempos, pares)
abajo.set_ylabel("par (N·m)")
abajo.set_xlabel("tiempo (s)")
for eje in (arriba, medio, abajo):
    eje.axvspan(INICIO_EMPUJON, FIN_EMPUJON, color="gray", alpha=0.25, linewidth=0)
    eje.grid(alpha=0.3)
figura.tight_layout()

if len(sys.argv) > 1:
    figura.savefig(sys.argv[1], dpi=120)
else:
    plt.show()
