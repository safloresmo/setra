"""Control de posición de un motor de corriente continua.

Un perfil trapezoidal da la referencia y un PID calcula el voltaje, ayudado
por un adelanto que ya pone el voltaje que el modelo dice que hace falta. A
media carrera entra un par de carga para ver cómo lo corrige la integral.

    python ejemplos/motor_dc.py              # abre la gráfica
    python ejemplos/motor_dc.py motor.png    # la guarda en un archivo
"""

import sys

import matplotlib

if len(sys.argv) > 1:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt

from setra import MotorDc, PerfilTrapezoidal, Pid

DT = 1e-3
VOLTAJE_MAX = 12.0
PAR_CARGA = 0.03
INICIO_CARGA = 0.2

motor = MotorDc(resistencia=2.0, inductancia=2e-3, kt=0.05, ke=0.05, inercia=2e-4, friccion=1e-4)
# Voltaje por unidad de velocidad y de aceleración, sin carga.
kv = motor.ke + motor.resistencia * motor.friccion / motor.kt
ka = motor.resistencia * motor.inercia / motor.kt

perfil = PerfilTrapezoidal(velocidad=100.0, aceleracion=400.0, desde=0.0, hasta=20.0)
pid = Pid(
    13.0,
    120.0,
    0.55,
    salida_min=-VOLTAJE_MAX,
    salida_max=VOLTAJE_MAX,
    tau_derivada=1e-3,
    derivada_del_error=True,
)

tiempos, referencias, posiciones, velocidades, voltajes = [], [], [], [], []
for n in range(round((perfil.duracion + 0.4) / DT)):
    t = n * DT
    referencia, velocidad = perfil.estado(t)
    aceleracion = (perfil.estado(t + DT)[1] - velocidad) / DT
    voltaje = pid.actualizar(referencia, motor.posicion, DT, adelanto=kv * velocidad + ka * aceleracion)

    tiempos.append(t)
    referencias.append(referencia)
    posiciones.append(motor.posicion)
    velocidades.append(motor.velocidad)
    voltajes.append(voltaje)

    # La planta se integra más fino que el lazo de control.
    motor.paso(voltaje, DT, par_carga=PAR_CARGA if t >= INICIO_CARGA else 0.0, subpasos=20)

errores = [r - p for r, p in zip(referencias, posiciones)]
print(f"error máximo de seguimiento: {max(map(abs, errores)):.4f} rad")
print(f"error final: {abs(errores[-1]):.5f} rad")

figura, (arriba, medio, abajo) = plt.subplots(3, 1, sharex=True, figsize=(8, 7))
arriba.plot(tiempos, referencias, "--", label="referencia")
arriba.plot(tiempos, posiciones, label="posición")
arriba.set_ylabel("rad")
arriba.legend()
medio.plot(tiempos, errores)
medio.set_ylabel("error (rad)")
abajo.plot(tiempos, voltajes)
abajo.set_ylabel("voltaje (V)")
abajo.set_xlabel("tiempo (s)")
for eje in (arriba, medio, abajo):
    eje.axvline(INICIO_CARGA, color="gray", linewidth=0.8, linestyle=":")
    eje.grid(alpha=0.3)
figura.tight_layout()

if len(sys.argv) > 1:
    figura.savefig(sys.argv[1], dpi=120)
else:
    plt.show()
