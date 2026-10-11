//! Control de posición de un motor de corriente continua: un perfil
//! trapezoidal da la referencia y un PID calcula el voltaje, ayudado por un
//! adelanto que ya pone el voltaje que el modelo dice que hace falta. Escribe
//! un CSV por la salida estándar.
//!
//! ```sh
//! cargo run --example motor_dc > motor.csv
//! ```

use setra::Real;
use setra::perfil::{Estado, Limites, PerfilTrapezoidal};
use setra::pid::Pid;
use setra::planta::MotorDc;

const DT: Real = 1e-3;
// La planta se integra más fino que el lazo de control.
const SUBPASOS: usize = 20;
const VOLTAJE_MAX: Real = 12.0;

fn main() {
    let mut motor = MotorDc {
        resistencia: 2.0,
        inductancia: 2e-3,
        kt: 0.05,
        ke: 0.05,
        inercia: 2e-4,
        friccion: 1e-4,
        corriente: 0.0,
        velocidad: 0.0,
        posicion: 0.0,
    };
    // Voltaje por unidad de velocidad y de aceleración, sin carga.
    let kv = motor.ke + motor.resistencia * motor.friccion / motor.kt;
    let ka = motor.resistencia * motor.inercia / motor.kt;

    let perfil = PerfilTrapezoidal::planear(
        Limites {
            velocidad: 100.0,
            aceleracion: 400.0,
        },
        Estado::default(),
        Estado {
            posicion: 20.0,
            velocidad: 0.0,
        },
    );
    let mut pid = Pid::nuevo(13.0, 120.0, 0.55)
        .con_limites(-VOLTAJE_MAX, VOLTAJE_MAX)
        .con_filtro_derivada(1e-3)
        .con_derivada_del_error();

    println!("t,referencia,posicion,velocidad,voltaje");
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
        println!(
            "{t:.3},{:.5},{:.5},{:.4},{voltaje:.4}",
            referencia.posicion, motor.posicion, motor.velocidad
        );
        for _ in 0..SUBPASOS {
            motor.paso(voltaje, 0.0, DT / SUBPASOS as Real);
        }
    }
}
