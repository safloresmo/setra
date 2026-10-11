import math

import pytest

from setra import (
    Bicuadratico,
    LimitadorPendiente,
    MotorDc,
    PasoBajo,
    Pendulo,
    PerfilTrapezoidal,
    Pid,
)

DT = 1e-3


def test_pid_proporcional():
    assert Pid(2.0).actualizar(1.0, 0.25, DT) == 1.5


def test_pid_limites_y_adelanto():
    pid = Pid(1.0, salida_min=-12.0, salida_max=12.0)
    assert pid.actualizar(1.0, 0.0, DT, adelanto=3.0) == 4.0
    assert pid.actualizar(100.0, 0.0, DT) == 12.0


def test_pid_integral_no_crece_saturada():
    pid = Pid(1.0, 50.0, salida_min=-1.0, salida_max=1.0)
    for _ in range(1000):
        pid.actualizar(10.0, 0.0, DT)
    assert pid.integral == 0.0


def test_pid_ganancias_modificables():
    pid = Pid(1.0)
    pid.kp = 3.0
    assert pid.kp == 3.0
    assert pid.actualizar(1.0, 0.0, DT) == 3.0


def test_pid_derivada_del_error():
    sobre_medida = Pid(0.0, 0.0, 1.0)
    sobre_error = Pid(0.0, 0.0, 1.0, derivada_del_error=True)
    for n in range(10):
        x = n * 0.5
        a = sobre_medida.actualizar(x, x, 0.5)
        b = sobre_error.actualizar(x, x, 0.5)
    assert (a, b) == (-1.0, 0.0)


def test_paso_bajo():
    filtro = PasoBajo(0.1)
    filtro.actualizar(0.0, 1e-4)
    for _ in range(1000):
        salida = filtro.actualizar(1.0, 1e-4)
    assert salida == pytest.approx(1 - math.exp(-1), abs=1e-3)
    assert PasoBajo.con_corte(10.0).actualizar(3.0, DT) == 3.0


def amplitud(filtro, frecuencia, muestreo=1000.0):
    suma = 0.0
    for n in range(4000):
        salida = filtro.actualizar(math.sin(2 * math.pi * frecuencia * n / muestreo))
        if n >= 2000:
            suma += salida * salida
    return math.sqrt(2 * suma / 2000)


def test_bicuadratico():
    assert amplitud(Bicuadratico.paso_bajo(50.0, 1000.0), 50.0) == pytest.approx(math.sqrt(0.5), abs=1e-3)
    assert amplitud(Bicuadratico.rechaza_banda(60.0, 1000.0, 5.0), 60.0) < 1e-3
    assert Bicuadratico([1.0, 0.0, 0.0], [0.0, 0.0]).actualizar(7.0) == 7.0


def test_limitador_de_pendiente():
    limitador = LimitadorPendiente(2.0)
    limitador.actualizar(0.0, 0.1)
    assert limitador.actualizar(10.0, 0.1) == pytest.approx(0.2)


def test_perfil_trapezoidal():
    perfil = PerfilTrapezoidal(2.0, 4.0, 0.0, 5.0)
    assert perfil.duracion == pytest.approx(3.0)
    assert perfil.estado(1.5)[1] == 2.0
    assert perfil.estado(10.0) == (5.0, 0.0)
    assert not perfil.terminado(2.9) and perfil.terminado(3.0)

    medio = perfil.estado(1.0)
    seguido = PerfilTrapezoidal(2.0, 4.0, medio[0], 8.0, velocidad_inicial=medio[1])
    assert seguido.estado(0.0) == medio


def test_motor_velocidad_de_regimen():
    motor = MotorDc(2.0, 2e-3, 0.05, 0.05, 2e-4, 1e-4)
    for _ in range(2000):
        motor.paso(12.0, DT, subpasos=20)
    esperada = motor.kt * 12.0 / (motor.resistencia * motor.friccion + motor.kt * motor.ke)
    assert motor.velocidad == pytest.approx(esperada, rel=1e-5)
    assert motor.posicion > 0.0 and motor.corriente > 0.0


def test_pendulo_conserva_la_energia():
    pendulo = Pendulo(0.2, 0.5, angulo=1.0)
    inicial = pendulo.energia
    for _ in range(10_000):
        pendulo.paso(0.0, DT)
    assert pendulo.energia == pytest.approx(inicial, rel=1e-6)
    assert pendulo.velocidad != 0.0 and abs(pendulo.angulo) <= 1.0 + 1e-6
