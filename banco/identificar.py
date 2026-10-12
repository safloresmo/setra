"""Banco de identificación de motores.

Aplica al motor una serie de escalones de voltaje, registra la posición del
encoder y ajusta el modelo con `setra.identificar_motor`. Al terminar escribe
los parámetros del motor, las ganancias de un PI de velocidad y dibuja la
medición junto a lo que predice el modelo.

    python banco/identificar.py --simulado
    python banco/identificar.py --puerto /dev/ttyUSB0 --voltaje-fuente 12 --cuentas 1320

Con `--puerto` habla con el programa de `banco/firmware/`. Necesita
`matplotlib` y, para el banco real, `pyserial`.
"""

import argparse
import csv
import math
import sys

from setra import ModeloMotor, MotorDc, identificar_motor

PERIODO = 1e-3  # segundos entre muestras


def opciones():
    o = argparse.ArgumentParser(description="Identifica un motor de corriente continua con encoder.")
    origen = o.add_mutually_exclusive_group(required=True)
    origen.add_argument("--puerto", help="puerto serie del banco, por ejemplo /dev/ttyUSB0")
    origen.add_argument("--simulado", action="store_true", help="usa un motor simulado en vez del banco")
    origen.add_argument("--csv", help="vuelve a analizar un ensayo guardado")
    o.add_argument("--voltaje-fuente", type=float, default=12.0, help="voltaje de la fuente del motor (V)")
    o.add_argument("--cuentas", type=float, default=1320.0, help="cuentas del encoder por vuelta, ya en cuadratura")
    o.add_argument("--niveles", type=float, nargs="+", default=[3.0, 6.0, 9.0, 12.0], help="voltajes de los escalones (V)")
    o.add_argument("--encendido", type=float, default=0.6, help="segundos que dura cada escalón")
    o.add_argument("--apagado", type=float, default=0.4, help="segundos de descanso tras cada escalón")
    o.add_argument("--resistencia", type=float, help="resistencia del devanado (ohmios), para calcular la inercia")
    o.add_argument("--tau-lazo", type=float, help="constante de tiempo deseada del lazo de velocidad (s)")
    o.add_argument("--guardar", help="guarda el ensayo en este CSV")
    o.add_argument("--grafica", help="guarda la gráfica en este archivo en vez de abrirla")
    return o.parse_args()


def programa(o):
    """Los tramos del ensayo: (milésimas del voltaje de la fuente, milisegundos)."""
    tramos = []
    for nivel in o.niveles:
        if abs(nivel) > o.voltaje_fuente:
            sys.exit(f"El nivel de {nivel} V pasa del voltaje de la fuente ({o.voltaje_fuente} V).")
        tramos.append((round(1000 * nivel / o.voltaje_fuente), round(1000 * o.encendido)))
        tramos.append((0, round(1000 * o.apagado)))
    return tramos


def ensayo_simulado(tramos, o):
    """Un motor con los datos de un motorreductor pequeño y un encoder ideal."""
    motor = MotorDc(2.0, 2e-3, 0.05, 0.05, 2e-4, 1e-4, friccion_seca=0.01)
    ciclos, cuentas = [], []
    for milesimas, duracion in tramos:
        for _ in range(duracion):
            ciclos.append(milesimas)
            cuentas.append(round(motor.posicion * o.cuentas / (2 * math.pi)))
            motor.paso(milesimas / 1000 * o.voltaje_fuente, PERIODO, subpasos=20)
    return ciclos, cuentas


def ensayo_real(tramos, o):
    import serial

    orden = f"ENSAYO {round(PERIODO * 1e6)} {len(tramos)} " + " ".join(f"{m} {d}" for m, d in tramos)
    duracion = sum(d for _, d in tramos) / 1000
    with serial.Serial(o.puerto, 115200, timeout=duracion + 10) as puerto:
        puerto.reset_input_buffer()
        puerto.write((orden + "\n").encode())
        print(f"Ensayo en marcha, {duracion:.1f} s...")
        ciclos, cuentas = [], []
        while True:
            linea = puerto.readline().decode(errors="replace").strip()
            if not linea:
                sys.exit("El banco dejó de responder.")
            if linea.startswith("ERROR"):
                sys.exit(f"El banco rechazó el ensayo: {linea}")
            if linea == "FIN":
                return ciclos, cuentas
            if "," in linea:
                milesimas, cuenta = linea.split(",")
                ciclos.append(int(milesimas))
                cuentas.append(int(cuenta))


def leer_csv(ruta):
    with open(ruta) as archivo:
        filas = list(csv.DictReader(archivo))
    return [float(f["voltaje"]) for f in filas], [float(f["posicion"]) for f in filas]


def simular(modelo, voltajes):
    """Posición y velocidad que predice el modelo para esos voltajes."""
    posicion, velocidad, posiciones, velocidades = 0.0, 0.0, [], []
    pasos = 10
    dt = PERIODO / pasos
    for v in voltajes:
        posiciones.append(posicion)
        velocidades.append(velocidad)
        for _ in range(pasos):
            friccion = modelo.voltaje_friccion * math.tanh(velocidad / 0.5)
            velocidad += dt * (modelo.ganancia * (v - friccion) - velocidad) / modelo.tau
            posicion += dt * velocidad
    return posiciones, velocidades


def velocidad_medida(posiciones, ventana=10):
    """Velocidad por diferencias, promediada para que el encoder no la llene de ruido."""
    n = len(posiciones)
    return [
        (posiciones[min(k + ventana, n - 1)] - posiciones[max(k - ventana, 0)])
        / ((min(k + ventana, n - 1) - max(k - ventana, 0)) * PERIODO)
        for k in range(n)
    ]


def main():
    o = opciones()

    if o.csv:
        voltajes, posiciones = leer_csv(o.csv)
    else:
        tramos = programa(o)
        ciclos, cuentas = ensayo_simulado(tramos, o) if o.simulado else ensayo_real(tramos, o)
        voltajes = [m / 1000 * o.voltaje_fuente for m in ciclos]
        posiciones = [(c - cuentas[0]) * 2 * math.pi / o.cuentas for c in cuentas]

    if o.guardar:
        with open(o.guardar, "w", newline="") as archivo:
            escritor = csv.writer(archivo)
            escritor.writerow(["tiempo", "voltaje", "posicion"])
            escritor.writerows((f"{k * PERIODO:.3f}", v, f"{p:.6f}") for k, (v, p) in enumerate(zip(voltajes, posiciones)))

    try:
        modelo = identificar_motor(voltajes, posiciones, PERIODO)
    except ValueError as error:
        sys.exit(f"No se pudo identificar el motor: {error}. ¿Giró el motor y cuenta el encoder?")

    previstas, velocidades_modelo = simular(modelo, voltajes)
    error = math.sqrt(sum((a - b) ** 2 for a, b in zip(posiciones, previstas)) / len(posiciones))
    recorrido = max(posiciones) - min(posiciones)

    print(f"Ganancia:              {modelo.ganancia:.2f} rad/s por V")
    print(f"Constante de tiempo:   {modelo.tau * 1000:.1f} ms")
    print(f"Voltaje de fricción:   {modelo.voltaje_friccion:.2f} V")
    print(f"Constante del motor:   {1 / modelo.ganancia:.4f} V·s/rad (aproximada)")
    if o.resistencia:
        print(f"Inercia:               {modelo.inercia(o.resistencia):.2e} kg·m²")
    tau_lazo = o.tau_lazo or modelo.tau / 3
    kp, ki = modelo.pi_velocidad(tau_lazo)
    print(f"PI de velocidad:       kp = {kp:.4f}, ki = {ki:.3f}  (lazo de {tau_lazo * 1000:.0f} ms)")
    print(f"Error del modelo:      {error:.2f} rad en un recorrido de {recorrido:.0f} rad ({100 * error / recorrido:.2f} %)")

    import matplotlib

    if o.grafica:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    tiempos = [k * PERIODO for k in range(len(voltajes))]
    figura, (arriba, abajo) = plt.subplots(2, 1, sharex=True, figsize=(8, 5.5))
    arriba.plot(tiempos, velocidad_medida(posiciones), label="medida")
    arriba.plot(tiempos, velocidades_modelo, "--", label="modelo")
    arriba.set_ylabel("velocidad (rad/s)")
    arriba.legend()
    abajo.plot(tiempos, voltajes, color="gray")
    abajo.set_ylabel("voltaje (V)")
    abajo.set_xlabel("tiempo (s)")
    for eje in (arriba, abajo):
        eje.grid(alpha=0.3)
    figura.tight_layout()
    if o.grafica:
        figura.savefig(o.grafica, dpi=120)
    else:
        plt.show()


if __name__ == "__main__":
    main()
