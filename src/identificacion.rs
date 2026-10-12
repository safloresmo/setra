//! Identificación: ajustar un modelo a datos medidos.

use crate::Real;

/// Modelo de velocidad de un motor de corriente continua:
///
/// ```text
/// tau·dω/dt + ω = ganancia·(V − voltaje_friccion·signo(ω))
/// ```
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ModeloMotor {
    /// Velocidad de régimen por voltio, en rad/s por V.
    pub ganancia: Real,
    /// Constante de tiempo mecánica, en segundos.
    pub tau: Real,
    /// Voltaje que se va en vencer la fricción seca, en V.
    pub voltaje_friccion: Real,
}

impl ModeloMotor {
    /// Ganancias `(kp, ki)` de un PI de velocidad que deja el lazo cerrado
    /// como un sistema de primer orden con constante de tiempo `tau_lazo`.
    pub fn pi_velocidad(&self, tau_lazo: Real) -> (Real, Real) {
        let ki = 1.0 / (self.ganancia * tau_lazo);
        (self.tau * ki, ki)
    }

    /// Inercia del rotor y la carga, en kg·m², a partir de la resistencia del
    /// devanado. Supone que la fricción viscosa es pequeña frente al freno
    /// eléctrico, que es lo habitual en motores pequeños.
    pub fn inercia(&self, resistencia: Real) -> Real {
        self.tau / (resistencia * self.ganancia * self.ganancia)
    }
}

/// Ajusta un [`ModeloMotor`] a un ensayo: el voltaje aplicado y la posición
/// del eje (rad), muestreados cada `dt` segundos.
///
/// Trabaja con la posición y no con la velocidad, porque un encoder da la
/// posición exacta salvo su resolución, mientras que la velocidad calculada
/// por diferencias sale muy ruidosa. Para separar la fricción seca de la
/// ganancia el ensayo debe tener más de un nivel de voltaje.
///
/// Devuelve `None` si hay pocas muestras o si los datos no determinan el
/// modelo, por ejemplo si el motor no se movió.
// Las conversiones a f64 sobran cuando `Real` ya lo es.
#[allow(clippy::unnecessary_cast)]
pub fn identificar_motor(voltaje: &[Real], posicion: &[Real], dt: Real) -> Option<ModeloMotor> {
    let n = voltaje.len().min(posicion.len());
    if n < 8 || dt <= 0.0 {
        return None;
    }
    let dt = dt as f64;
    let inicio = posicion[0] as f64;

    // Integrando dos veces la ecuación del modelo queda una relación lineal
    // entre la posición y cuatro señales que se calculan de los datos:
    //   θ − θ0 = p0·∬V + p1·∫(θ − θ0) + p2·∬signo(ω) + p3·t
    // Los coeficientes se ajustan por mínimos cuadrados. Las sumas van en
    // f64: crecen mucho y en f32 se perdería precisión.
    let (mut int_v, mut int2_v) = (0.0, 0.0);
    let (mut int_s, mut int2_s) = (0.0, 0.0);
    let mut int_pos = 0.0;
    let mut ata = [[0.0f64; 4]; 4];
    let mut atb = [0.0f64; 4];

    for k in 1..n {
        let previa = posicion[k - 1] as f64 - inicio;
        let actual = posicion[k] as f64 - inicio;
        int_pos += (previa + actual) / 2.0 * dt;

        // El voltaje se mantiene durante todo el intervalo.
        let v = voltaje[k - 1] as f64;
        int2_v += (int_v + v * dt / 2.0) * dt;
        int_v += v * dt;

        let avance = actual - previa;
        let signo = if avance > 0.0 {
            1.0
        } else if avance < 0.0 {
            -1.0
        } else {
            0.0
        };
        int2_s += (int_s + signo * dt / 2.0) * dt;
        int_s += signo * dt;

        let fila = [int2_v, int_pos, int2_s, k as f64 * dt];
        for i in 0..4 {
            atb[i] += fila[i] * actual;
            for j in 0..4 {
                ata[i][j] += fila[i] * fila[j];
            }
        }
    }

    // Con un solo nivel de voltaje la fricción no se distingue de la
    // ganancia: se repite el ajuste sin ella.
    let p = resolver(ata, atb, [true; 4])
        .filter(|p| modelo(p).is_some())
        .or_else(|| resolver(ata, atb, [true, true, false, true]))?;
    modelo(&p)
}

#[allow(clippy::unnecessary_cast)]
fn modelo(p: &[f64; 4]) -> Option<ModeloMotor> {
    let tau = -1.0 / p[1];
    let ganancia = p[0] * tau;
    let voltaje_friccion = -p[2] / p[0];
    let valido = tau.is_finite()
        && tau > 0.0
        && ganancia.is_finite()
        && ganancia != 0.0
        && voltaje_friccion.is_finite()
        && voltaje_friccion >= 0.0;
    valido.then_some(ModeloMotor {
        ganancia: ganancia as Real,
        tau: tau as Real,
        voltaje_friccion: voltaje_friccion as Real,
    })
}

// Resuelve el sistema por eliminación de Gauss con pivoteo, usando solo las
// incógnitas marcadas; las demás quedan en cero.
#[allow(clippy::needless_range_loop)]
fn resolver(mut a: [[f64; 4]; 4], mut b: [f64; 4], usar: [bool; 4]) -> Option<[f64; 4]> {
    for i in 0..4 {
        if !usar[i] {
            for j in 0..4 {
                a[i][j] = 0.0;
                a[j][i] = 0.0;
            }
            a[i][i] = 1.0;
            b[i] = 0.0;
        }
    }

    for col in 0..4 {
        let mut pivote = col;
        for fila in col + 1..4 {
            if a[fila][col].abs() > a[pivote][col].abs() {
                pivote = fila;
            }
        }
        let escala = a[col].iter().fold(0.0f64, |m, v| m.max(v.abs()));
        if a[pivote][col].abs() <= 1e-12 * escala.max(f64::MIN_POSITIVE) {
            return None;
        }
        a.swap(col, pivote);
        b.swap(col, pivote);
        for fila in col + 1..4 {
            let factor = a[fila][col] / a[col][col];
            for j in col..4 {
                a[fila][j] -= factor * a[col][j];
            }
            b[fila] -= factor * b[col];
        }
    }

    let mut x = [0.0f64; 4];
    for i in (0..4).rev() {
        let mut suma = b[i];
        for j in i + 1..4 {
            suma -= a[i][j] * x[j];
        }
        x[i] = suma / a[i][i];
    }
    Some(x)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::planta::MotorDc;

    const DT: Real = 1e-3;
    const PI: Real = core::f64::consts::PI as Real;

    fn motor(friccion_seca: Real) -> MotorDc {
        MotorDc {
            resistencia: 2.0,
            inductancia: 2e-3,
            kt: 0.05,
            ke: 0.05,
            inercia: 2e-4,
            friccion: 1e-4,
            friccion_seca,
            corriente: 0.0,
            velocidad: 0.0,
            posicion: 0.0,
        }
    }

    // Aplica cada voltaje durante 0.6 s y deja el motor suelto 0.4 s. La
    // posición se redondea a la resolución de un encoder.
    fn ensayar(
        mut m: MotorDc,
        niveles: &[Real],
        cuentas_por_vuelta: Real,
    ) -> (Vec<Real>, Vec<Real>) {
        let paso_encoder = 2.0 * PI / cuentas_por_vuelta;
        let (mut voltajes, mut posiciones) = (Vec::new(), Vec::new());
        for &nivel in niveles {
            for n in 0..1000 {
                let v = if n < 600 { nivel } else { 0.0 };
                voltajes.push(v);
                posiciones.push((m.posicion / paso_encoder).round() * paso_encoder);
                for _ in 0..20 {
                    m.paso(v, 0.0, DT / 20.0);
                }
            }
        }
        (voltajes, posiciones)
    }

    fn cerca(valor: Real, esperado: Real, tolerancia: Real) -> bool {
        (valor - esperado).abs() <= tolerancia * esperado.abs()
    }

    #[test]
    fn recupera_el_motor_sin_friccion_seca() {
        let m = motor(0.0);
        let (v, p) = ensayar(m, &[3.0, 6.0, 9.0, 12.0], 2000.0);
        let modelo = identificar_motor(&v, &p, DT).unwrap();

        let ganancia = m.kt / (m.resistencia * m.friccion + m.kt * m.ke);
        let tau = m.inercia * m.resistencia / (m.resistencia * m.friccion + m.kt * m.ke);
        assert!(cerca(modelo.ganancia, ganancia, 0.02), "{modelo:?}");
        assert!(cerca(modelo.tau, tau, 0.05), "{modelo:?}");
        assert!(modelo.voltaje_friccion < 0.05, "{modelo:?}");
    }

    #[test]
    fn separa_la_friccion_seca() {
        let m = motor(0.01);
        let (v, p) = ensayar(m, &[3.0, 6.0, 9.0, 12.0], 2000.0);
        let modelo = identificar_motor(&v, &p, DT).unwrap();

        let ganancia = m.kt / (m.resistencia * m.friccion + m.kt * m.ke);
        let friccion = m.friccion_seca * m.resistencia / m.kt;
        assert!(cerca(modelo.ganancia, ganancia, 0.03), "{modelo:?}");
        assert!(cerca(modelo.voltaje_friccion, friccion, 0.15), "{modelo:?}");
    }

    #[test]
    fn aguanta_un_encoder_de_poca_resolucion() {
        let m = motor(0.0);
        let (v, p) = ensayar(m, &[4.0, 8.0, 12.0], 44.0);
        let modelo = identificar_motor(&v, &p, DT).unwrap();
        let ganancia = m.kt / (m.resistencia * m.friccion + m.kt * m.ke);
        assert!(cerca(modelo.ganancia, ganancia, 0.03), "{modelo:?}");
    }

    #[test]
    fn un_solo_nivel_se_ajusta_sin_friccion() {
        let m = motor(0.0);
        let (v, p) = ensayar(m, &[6.0], 2000.0);
        let modelo = identificar_motor(&v[..600], &p[..600], DT).unwrap();
        let ganancia = m.kt / (m.resistencia * m.friccion + m.kt * m.ke);
        assert!(cerca(modelo.ganancia, ganancia, 0.03), "{modelo:?}");
        assert_eq!(modelo.voltaje_friccion, 0.0);
    }

    #[test]
    fn sin_movimiento_no_hay_modelo() {
        let voltajes = [0.0; 100];
        let posiciones = [0.0; 100];
        assert!(identificar_motor(&voltajes, &posiciones, DT).is_none());
        assert!(identificar_motor(&voltajes[..4], &posiciones[..4], DT).is_none());
    }

    #[test]
    fn el_pi_sugerido_da_el_lazo_pedido() {
        // Con ese PI el lazo cerrado es de primer orden: a un tau_lazo llega al 63 %.
        let modelo = ModeloMotor {
            ganancia: 18.5,
            tau: 0.15,
            voltaje_friccion: 0.0,
        };
        let (kp, ki) = modelo.pi_velocidad(0.05);
        let (mut w, mut integral): (Real, Real) = (0.0, 0.0);
        let dt = 1e-4;
        for _ in 0..500 {
            let error = 1.0 - w;
            integral += ki * error * dt;
            let u = kp * error + integral;
            w += dt * (modelo.ganancia * u - w) / modelo.tau;
        }
        assert!((w - 0.632).abs() < 0.01, "w = {w}");
    }
}
