//! Modelos de plantas para probar controladores sin hardware.

use crate::{Mate, Real};

/// Un paso de Runge-Kutta de orden 4 para `dx/dt = f(x)`. La entrada de la
/// planta va capturada en `f` y se supone constante durante el paso.
pub fn rk4<const N: usize>(
    f: impl Fn(&[Real; N]) -> [Real; N],
    x: &[Real; N],
    dt: Real,
) -> [Real; N] {
    let avanzar = |k: &[Real; N], h: Real| core::array::from_fn(|i| x[i] + h * k[i]);

    let k1 = f(x);
    let k2 = f(&avanzar(&k1, dt / 2.0));
    let k3 = f(&avanzar(&k2, dt / 2.0));
    let k4 = f(&avanzar(&k3, dt));

    core::array::from_fn(|i| x[i] + dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]))
}

/// Motor de corriente continua de imanes permanentes, en unidades del SI.
#[derive(Clone, Copy, Debug)]
pub struct MotorDc {
    /// Resistencia del devanado, en ohmios.
    pub resistencia: Real,
    /// Inductancia del devanado, en henrios.
    pub inductancia: Real,
    /// Constante de par, en N·m/A.
    pub kt: Real,
    /// Constante de fuerza contraelectromotriz, en V·s/rad.
    pub ke: Real,
    /// Inercia del rotor y de la carga, en kg·m².
    pub inercia: Real,
    /// Fricción viscosa, en N·m·s/rad.
    pub friccion: Real,

    pub corriente: Real,
    pub velocidad: Real,
    pub posicion: Real,
}

impl MotorDc {
    /// Avanza `dt` segundos con el voltaje y el par de carga dados. La parte
    /// eléctrica es rápida: `dt` debe quedar por debajo de
    /// `inductancia / resistencia` o la integración se vuelve inestable.
    pub fn paso(&mut self, voltaje: Real, par_carga: Real, dt: Real) {
        let m = *self;
        let derivada = |x: &[Real; 3]| {
            let [corriente, velocidad, _] = *x;
            [
                (voltaje - m.resistencia * corriente - m.ke * velocidad) / m.inductancia,
                (m.kt * corriente - m.friccion * velocidad - par_carga) / m.inercia,
                velocidad,
            ]
        };
        let x = rk4(
            derivada,
            &[self.corriente, self.velocidad, self.posicion],
            dt,
        );
        [self.corriente, self.velocidad, self.posicion] = x;
    }
}

/// Péndulo simple con par aplicado en el eje. El ángulo se mide desde la
/// vertical hacia abajo, así que `PI` es la posición invertida.
#[derive(Clone, Copy, Debug)]
pub struct Pendulo {
    /// Masa concentrada en el extremo, en kg.
    pub masa: Real,
    /// Distancia del eje a la masa, en metros.
    pub longitud: Real,
    /// Fricción viscosa en el eje, en N·m·s/rad.
    pub friccion: Real,
    /// Aceleración de la gravedad, en m/s².
    pub gravedad: Real,

    pub angulo: Real,
    pub velocidad: Real,
}

impl Pendulo {
    pub fn paso(&mut self, par: Real, dt: Real) {
        let m = *self;
        let inercia = m.masa * m.longitud * m.longitud;
        let derivada = |x: &[Real; 2]| {
            let [angulo, velocidad] = *x;
            [
                velocidad,
                (par - m.friccion * velocidad) / inercia
                    - m.gravedad / m.longitud * Mate::sin(angulo),
            ]
        };
        [self.angulo, self.velocidad] = rk4(derivada, &[self.angulo, self.velocidad], dt);
    }

    /// Energía mecánica, con el cero en el punto más bajo.
    pub fn energia(&self) -> Real {
        let inercia = self.masa * self.longitud * self.longitud;
        inercia * self.velocidad * self.velocidad / 2.0
            + self.masa * self.gravedad * self.longitud * (1.0 - Mate::cos(self.angulo))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::PI;

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

    fn pendulo() -> Pendulo {
        Pendulo {
            masa: 0.2,
            longitud: 0.5,
            friccion: 0.0,
            gravedad: 9.81,
            angulo: 0.0,
            velocidad: 0.0,
        }
    }

    #[test]
    fn rk4_sigue_una_exponencial() {
        let mut x = [1.0];
        for _ in 0..100 {
            x = rk4(|x| [-x[0]], &x, 0.01);
        }
        assert!((x[0] - Mate::exp(-1.0)).abs() < 1e-5);
    }

    #[test]
    fn el_motor_llega_a_su_velocidad_de_regimen() {
        let mut m = motor();
        for _ in 0..40_000 {
            m.paso(12.0, 0.0, 5e-5);
        }
        let esperada = m.kt * 12.0 / (m.resistencia * m.friccion + m.kt * m.ke);
        assert!(
            (m.velocidad - esperada).abs() / esperada < 1e-3,
            "w = {}",
            m.velocidad
        );
        // En régimen el par del motor solo vence la fricción.
        assert!((m.kt * m.corriente - m.friccion * m.velocidad).abs() < 1e-4);
    }

    #[test]
    fn el_par_de_carga_frena_el_motor() {
        let (mut libre, mut cargado) = (motor(), motor());
        for _ in 0..40_000 {
            libre.paso(12.0, 0.0, 5e-5);
            cargado.paso(12.0, 0.1, 5e-5);
        }
        assert!(cargado.velocidad < libre.velocidad - 50.0);
    }

    #[test]
    fn el_pendulo_sin_friccion_conserva_la_energia() {
        let mut p = pendulo();
        p.angulo = 1.0;
        let inicial = p.energia();
        for _ in 0..10_000 {
            p.paso(0.0, 1e-3);
        }
        assert!((p.energia() - inicial).abs() / inicial < 1e-3);
    }

    #[test]
    fn el_periodo_de_oscilaciones_pequenas() {
        let mut p = pendulo();
        p.angulo = 0.05;
        // Medio periodo después está en el extremo contrario.
        let medio_periodo = PI * Mate::sqrt(p.longitud / p.gravedad);
        let pasos = 5000;
        for _ in 0..pasos {
            p.paso(0.0, medio_periodo / pasos as Real);
        }
        assert!((p.angulo + 0.05).abs() < 1e-4, "angulo = {}", p.angulo);
    }
}
