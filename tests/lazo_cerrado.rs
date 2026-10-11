//! Lazos completos sobre el modelo del motor.

use setra::Real;
use setra::estados::{Observador, discretizar, ganancia_observador, lqr};
use setra::matriz::{Matriz, Vector};
use setra::perfil::{Estado, Limites, PerfilTrapezoidal};
use setra::pid::Pid;
use setra::planta::{MotorDc, Pendulo};

const PI: Real = core::f64::consts::PI as Real;

const DT: Real = 1e-3;
const SUBPASOS: usize = 20;
const VOLTAJE_MAX: Real = 12.0;

fn motor() -> MotorDc {
    MotorDc {
        resistencia: 2.0,
        inductancia: 2e-3,
        kt: 0.05,
        ke: 0.05,
        inercia: 2e-4,
        friccion: 1e-4,
        corriente: 0.0,
        velocidad: 0.0,
        posicion: 0.0,
    }
}

fn avanzar(motor: &mut MotorDc, voltaje: Real, par_carga: Real) {
    for _ in 0..SUBPASOS {
        motor.paso(voltaje, par_carga, DT / SUBPASOS as Real);
    }
}

#[test]
fn la_posicion_sigue_el_perfil() {
    let mut motor = motor();
    let kv = motor.ke + motor.resistencia * motor.friccion / motor.kt;
    let ka = motor.resistencia * motor.inercia / motor.kt;

    let meta = Estado {
        posicion: 20.0,
        velocidad: 0.0,
    };
    let limites = Limites {
        velocidad: 100.0,
        aceleracion: 400.0,
    };
    let perfil = PerfilTrapezoidal::planear(limites, Estado::default(), meta);
    let mut pid = Pid::nuevo(13.0, 120.0, 0.55)
        .con_limites(-VOLTAJE_MAX, VOLTAJE_MAX)
        .con_filtro_derivada(1e-3)
        .con_derivada_del_error();

    let mut peor_error: Real = 0.0;
    let pasos = ((perfil.duracion() + 0.3) / DT) as usize;
    for n in 0..pasos {
        let t = n as Real * DT;
        let referencia = perfil.estado(t);
        let aceleracion = (perfil.estado(t + DT).velocidad - referencia.velocidad) / DT;
        let voltaje = pid.actualizar_con_adelanto(
            referencia.posicion,
            motor.posicion,
            kv * referencia.velocidad + ka * aceleracion,
            DT,
        );
        assert!(voltaje.abs() <= VOLTAJE_MAX);
        peor_error = peor_error.max((referencia.posicion - motor.posicion).abs());
        avanzar(&mut motor, voltaje, 0.0);
    }

    assert!(peor_error < 0.05, "error de seguimiento = {peor_error}");
    assert!(
        (motor.posicion - meta.posicion).abs() < 5e-3,
        "posicion = {}",
        motor.posicion
    );
    assert!(
        motor.velocidad.abs() < 0.1,
        "velocidad = {}",
        motor.velocidad
    );
}

#[test]
fn la_velocidad_aguanta_un_par_de_carga() {
    const REFERENCIA: Real = 150.0;

    let mut motor = motor();
    let mut pid = Pid::nuevo(0.2, 4.0, 0.0).con_limites(-VOLTAJE_MAX, VOLTAJE_MAX);

    let mut maxima: Real = 0.0;
    let mut minima_con_carga = REFERENCIA;
    for n in 0..2000 {
        // A mitad de la prueba entra una carga de 0.05 N·m.
        let par_carga = if n >= 1000 { 0.05 } else { 0.0 };
        let voltaje = pid.actualizar(REFERENCIA, motor.velocidad, DT);
        assert!(voltaje.abs() <= VOLTAJE_MAX);
        avanzar(&mut motor, voltaje, par_carga);

        maxima = maxima.max(motor.velocidad);
        if n >= 1000 {
            minima_con_carga = minima_con_carga.min(motor.velocidad);
        }
        if n == 999 {
            assert!(
                (motor.velocidad - REFERENCIA).abs() < 0.5,
                "w = {}",
                motor.velocidad
            );
        }
    }

    // Arranca saturado: sin protección la integral acumulada lo pasaría de largo.
    assert!(maxima < REFERENCIA * 1.02, "maxima = {maxima}");
    assert!(
        minima_con_carga > REFERENCIA * 0.9,
        "minima = {minima_con_carga}"
    );
    assert!(
        (motor.velocidad - REFERENCIA).abs() < 0.5,
        "w = {}",
        motor.velocidad
    );
}

#[test]
fn el_pendulo_invertido_se_sostiene_midiendo_solo_el_angulo() {
    const PERIODO: Real = 2e-3;
    const PAR_MAX: Real = 0.5;

    let mut pendulo = Pendulo {
        masa: 0.2,
        longitud: 0.5,
        friccion: 1e-3,
        gravedad: 9.81,
        angulo: PI + 0.2,
        velocidad: 0.0,
    };

    // Modelo lineal alrededor de la vertical hacia arriba; el estado es la
    // desviación del ángulo y su velocidad.
    let inercia = pendulo.masa * pendulo.longitud * pendulo.longitud;
    let (a, b) = discretizar(
        &Matriz([
            [0.0, 1.0],
            [
                pendulo.gravedad / pendulo.longitud,
                -pendulo.friccion / inercia,
            ],
        ]),
        &Matriz([[0.0], [1.0 / inercia]]),
        PERIODO,
    );
    let c = Matriz([[1.0, 0.0]]);
    let k = lqr(&a, &b, &Matriz::diagonal([10.0, 0.1]), &Matriz([[1.0]])).unwrap();
    let l =
        ganancia_observador(&a, &c, &Matriz::diagonal([1.0, 100.0]), &Matriz([[0.01]])).unwrap();
    let mut observador = Observador {
        a,
        b,
        c,
        l,
        x: Vector::ceros(),
    };

    for _ in 0..1500 {
        let par = (-(k * observador.x).0[0][0]).clamp(-PAR_MAX, PAR_MAX);
        observador.actualizar(
            &Vector::columna([par]),
            &Vector::columna([pendulo.angulo - PI]),
        );
        for _ in 0..4 {
            pendulo.paso(par, PERIODO / 4.0);
        }
    }

    assert!(
        (pendulo.angulo - PI).abs() < 5e-3,
        "angulo = {}",
        pendulo.angulo
    );
    assert!(
        pendulo.velocidad.abs() < 5e-2,
        "velocidad = {}",
        pendulo.velocidad
    );
    let estimada = observador.x.0[1][0];
    assert!(
        (estimada - pendulo.velocidad).abs() < 5e-2,
        "estimada = {estimada}"
    );
}
