//! Lazos completos sobre el modelo del motor.

use setra::Real;
use setra::perfil::{Estado, Limites, PerfilTrapezoidal};
use setra::pid::Pid;
use setra::planta::MotorDc;

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
