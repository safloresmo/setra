"""Péndulo de Furuta sostenido hacia arriba, con las limitaciones del montaje.

Sirve para dimensionar el montaje antes de construirlo: el motor se manda por
voltaje y satura, los ángulos llegan con la resolución de los encoders y las
velocidades las reconstruye un observador.

Los datos del motor, de los encoders y de las piezas son supuestos razonables
para un motorreductor de 12 V. Hay que cambiarlos por los medidos.

    python ejemplos/furuta.py               # abre la gráfica
    python ejemplos/furuta.py furuta.png    # la guarda en un archivo
"""

import math
import sys

import matplotlib

if len(sys.argv) > 1:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt

from setra import Furuta, Observador, discretizar, ganancia_observador, lqr

DT = 2e-3

# Motor, visto en el eje de salida de la reductora.
VOLTAJE_MAX = 12.0
RESISTENCIA = 4.0  # ohmios
KT = 0.25  # N·m/A; con el mismo valor en V·s/rad para la fuerza contraelectromotriz
FRICCION_MOTOR = 2e-3  # N·m·s/rad

# Cuentas por vuelta de cada encoder, ya en cuadratura.
CUENTAS_BRAZO = 924
CUENTAS_PENDULO = 2400

# Piezas: brazo de 12 cm y péndulo de barra de 20 cm y 20 g.
INERCIA_BRAZO = 2e-3
LONGITUD_BRAZO = 0.12
MASA_PENDULO = 0.02
LARGO_PENDULO = 0.20
FRICCION_PENDULO = 1e-4

DESVIACION_INICIAL = math.radians(8.0)
INICIO_EMPUJON, FIN_EMPUJON, PAR_EMPUJON = 2.0, 2.05, 0.05


def planta(friccion_brazo, **estado):
    return Furuta(
        INERCIA_BRAZO,
        LONGITUD_BRAZO,
        MASA_PENDULO,
        centro_pendulo=LARGO_PENDULO / 2,
        inercia_pendulo=MASA_PENDULO * LARGO_PENDULO**2 / 3,
        friccion_brazo=friccion_brazo,
        friccion_pendulo=FRICCION_PENDULO,
        **estado,
    )


def medir(angulo, cuentas):
    paso = 2 * math.pi / cuentas
    return round(angulo / paso) * paso


# Mandar el motor por voltaje equivale a mandarlo por par con una ganancia
# KT/R y con más fricción, la que pone la fuerza contraelectromotriz.
a, b = planta(FRICCION_MOTOR + KT * KT / RESISTENCIA).lineal_arriba()
b = [[fila[0] * KT / RESISTENCIA] for fila in b]
a, b = discretizar(a, b, DT)
c = [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0]]

(k,) = lqr(a, b, [[5, 0, 0, 0], [0, 50, 0, 0], [0, 0, 0.1, 0], [0, 0, 0, 0.5]], [[1.0]])
observador = Observador(
    a,
    b,
    c,
    ganancia_observador(a, c, [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 100, 0], [0, 0, 0, 100]], [[0.01, 0], [0, 0.01]]),
)

furuta = planta(FRICCION_MOTOR, angulo_pendulo=math.pi + DESVIACION_INICIAL)
tiempos, brazos, pendulos, voltajes = [], [], [], []
for n in range(round(4.0 / DT)):
    t = n * DT
    voltaje = -sum(ganancia * x for ganancia, x in zip(k, observador.x))
    voltaje = max(-VOLTAJE_MAX, min(VOLTAJE_MAX, voltaje))
    medida = [
        medir(furuta.angulo_brazo, CUENTAS_BRAZO),
        medir(furuta.angulo_pendulo - math.pi, CUENTAS_PENDULO),
    ]
    observador.actualizar(voltaje, medida)

    tiempos.append(t)
    brazos.append(math.degrees(furuta.angulo_brazo))
    pendulos.append(math.degrees(furuta.angulo_pendulo - math.pi))
    voltajes.append(voltaje)

    par = KT * (voltaje - KT * furuta.velocidad_brazo) / RESISTENCIA
    empujon = PAR_EMPUJON if INICIO_EMPUJON <= t < FIN_EMPUJON else 0.0
    furuta.paso(par + empujon, DT, subpasos=4)

final = round(3.5 / DT)
print("ganancias (V/rad, V/rad, V·s/rad, V·s/rad): " + ", ".join(f"{g:.2f}" for g in k))
print(f"recorrido máximo del brazo: {max(map(abs, brazos)):.1f}°")
print(f"voltaje máximo: {max(map(abs, voltajes)):.1f} V de {VOLTAJE_MAX:.0f}")
print(f"temblor final del péndulo: ±{max(map(abs, pendulos[final:])):.2f}°")
print(f"temblor final del brazo: ±{max(map(abs, brazos[final:])):.2f}°")

figura, (arriba, medio, abajo) = plt.subplots(3, 1, sharex=True, figsize=(8, 7))
arriba.plot(tiempos, pendulos)
arriba.set_ylabel("péndulo (°)")
medio.plot(tiempos, brazos)
medio.set_ylabel("brazo (°)")
abajo.plot(tiempos, voltajes)
abajo.set_ylabel("voltaje (V)")
abajo.set_xlabel("tiempo (s)")
for eje in (arriba, medio, abajo):
    eje.axvspan(INICIO_EMPUJON, FIN_EMPUJON, color="gray", alpha=0.25, linewidth=0)
    eje.grid(alpha=0.3)
figura.tight_layout()

if len(sys.argv) > 1:
    figura.savefig(sys.argv[1], dpi=120)
else:
    plt.show()
