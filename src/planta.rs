//! Modelos de plantas para probar controladores sin hardware.

use crate::matriz::Matriz;
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

// La fricción seca cambia de signo con el giro. Para poder integrarla se
// suaviza el cambio alrededor de cero, en un margen de esta velocidad (rad/s).
const VELOCIDAD_FRICCION_SECA: Real = 0.5;

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
    /// Fricción seca, en N·m: el par que hay que vencer para que gire.
    pub friccion_seca: Real,

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
                (m.kt * corriente
                    - m.friccion * velocidad
                    - m.friccion_seca * Mate::tanh(velocidad / VELOCIDAD_FRICCION_SECA)
                    - par_carga)
                    / m.inercia,
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

/// Péndulo invertido rotatorio, o de Furuta: un motor gira un brazo horizontal
/// que lleva en la punta un péndulo libre. El ángulo del péndulo se mide desde
/// la vertical hacia abajo, así que `PI` es la posición invertida.
///
/// El péndulo se trata como una barra delgada: solo cuenta su inercia respecto
/// al eje en el que gira.
#[derive(Clone, Copy, Debug)]
pub struct Furuta {
    /// Inercia del brazo y del rotor respecto al eje del motor, en kg·m².
    pub inercia_brazo: Real,
    /// Distancia del eje del motor al eje del péndulo, en metros.
    pub longitud_brazo: Real,
    /// Masa del péndulo, en kg.
    pub masa_pendulo: Real,
    /// Distancia del eje del péndulo a su centro de masa, en metros.
    pub centro_pendulo: Real,
    /// Inercia del péndulo respecto a su eje, en kg·m².
    pub inercia_pendulo: Real,
    /// Fricción viscosa en el eje del motor, en N·m·s/rad.
    pub friccion_brazo: Real,
    /// Fricción viscosa en el eje del péndulo, en N·m·s/rad.
    pub friccion_pendulo: Real,
    /// Aceleración de la gravedad, en m/s².
    pub gravedad: Real,

    pub angulo_brazo: Real,
    pub angulo_pendulo: Real,
    pub velocidad_brazo: Real,
    pub velocidad_pendulo: Real,
}

impl Furuta {
    /// Avanza `dt` segundos con el par dado en el eje del motor.
    pub fn paso(&mut self, par: Real, dt: Real) {
        let m = *self;
        let x = [
            self.angulo_brazo,
            self.angulo_pendulo,
            self.velocidad_brazo,
            self.velocidad_pendulo,
        ];
        let x = rk4(|x| m.derivada(x, par), &x, dt);
        [
            self.angulo_brazo,
            self.angulo_pendulo,
            self.velocidad_brazo,
            self.velocidad_pendulo,
        ] = x;
    }

    // Acoplamiento entre el brazo y el péndulo.
    fn acoplamiento(&self) -> Real {
        self.masa_pendulo * self.longitud_brazo * self.centro_pendulo
    }

    fn peso(&self) -> Real {
        self.masa_pendulo * self.gravedad * self.centro_pendulo
    }

    fn derivada(&self, x: &[Real; 4], par: Real) -> [Real; 4] {
        let [_, angulo, v_brazo, v_pendulo] = *x;
        let (seno, coseno) = (Mate::sin(angulo), Mate::cos(angulo));
        let jp = self.inercia_pendulo;

        // Matriz de inercia, simétrica de 2×2.
        let m11 = self.inercia_brazo
            + self.masa_pendulo * self.longitud_brazo * self.longitud_brazo
            + jp * seno * seno;
        let m12 = self.acoplamiento() * coseno;
        let m22 = jp;

        let f1 =
            par - self.friccion_brazo * v_brazo - 2.0 * jp * seno * coseno * v_brazo * v_pendulo
                + self.acoplamiento() * seno * v_pendulo * v_pendulo;
        let f2 = jp * seno * coseno * v_brazo * v_brazo
            - self.peso() * seno
            - self.friccion_pendulo * v_pendulo;

        let determinante = m11 * m22 - m12 * m12;
        [
            v_brazo,
            v_pendulo,
            (m22 * f1 - m12 * f2) / determinante,
            (m11 * f2 - m12 * f1) / determinante,
        ]
    }

    /// Energía mecánica, con el cero en el péndulo colgando y en reposo.
    pub fn energia(&self) -> Real {
        let (seno, coseno) = (
            Mate::sin(self.angulo_pendulo),
            Mate::cos(self.angulo_pendulo),
        );
        let m11 = self.inercia_brazo
            + self.masa_pendulo * self.longitud_brazo * self.longitud_brazo
            + self.inercia_pendulo * seno * seno;
        m11 * self.velocidad_brazo * self.velocidad_brazo / 2.0
            + self.acoplamiento() * coseno * self.velocidad_brazo * self.velocidad_pendulo
            + self.inercia_pendulo * self.velocidad_pendulo * self.velocidad_pendulo / 2.0
            + self.peso() * (1.0 - coseno)
    }

    /// Modelo lineal continuo `dx/dt = A·x + B·u` alrededor del péndulo hacia
    /// arriba. El estado es el ángulo del brazo, la desviación del péndulo
    /// respecto a la vertical y sus dos velocidades; la entrada es el par.
    pub fn lineal_arriba(&self) -> (Matriz<4, 4>, Matriz<4, 1>) {
        let j0 = self.inercia_brazo + self.masa_pendulo * self.longitud_brazo * self.longitud_brazo;
        let jp = self.inercia_pendulo;
        let c = self.acoplamiento();
        let d = j0 * jp - c * c;
        let (br, bp, peso) = (self.friccion_brazo, self.friccion_pendulo, self.peso());

        let a = Matriz([
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
            [0.0, c * peso / d, -jp * br / d, -c * bp / d],
            [0.0, j0 * peso / d, -c * br / d, -j0 * bp / d],
        ]);
        (a, Matriz([[0.0], [0.0], [jp / d], [c / d]]))
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
            friccion_seca: 0.0,
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

    fn furuta() -> Furuta {
        Furuta {
            inercia_brazo: 2e-3,
            longitud_brazo: 0.12,
            masa_pendulo: 0.02,
            centro_pendulo: 0.1,
            inercia_pendulo: 2.67e-4,
            friccion_brazo: 0.0,
            friccion_pendulo: 0.0,
            gravedad: 9.81,
            angulo_brazo: 0.0,
            angulo_pendulo: 0.0,
            velocidad_brazo: 0.0,
            velocidad_pendulo: 0.0,
        }
    }

    #[test]
    fn el_furuta_sin_friccion_conserva_la_energia() {
        let mut f = furuta();
        f.angulo_pendulo = 2.0;
        f.velocidad_brazo = 3.0;
        let inicial = f.energia();
        for _ in 0..10_000 {
            f.paso(0.0, 5e-4);
        }
        assert!((f.energia() - inicial).abs() / inicial < 1e-3);
        // El péndulo y el brazo se han intercambiado energía por el camino.
        assert!(
            (f.velocidad_brazo - 3.0).abs() > 0.05,
            "w = {}",
            f.velocidad_brazo
        );
    }

    #[test]
    fn el_furuta_colgando_oscila_como_un_pendulo() {
        // Con el brazo muy pesado queda un péndulo físico corriente.
        let mut f = furuta();
        f.inercia_brazo = 1e3;
        f.angulo_pendulo = 0.05;
        let peso = f.masa_pendulo * f.gravedad * f.centro_pendulo;
        let medio_periodo = PI * Mate::sqrt(f.inercia_pendulo / peso);
        let pasos = 5000;
        for _ in 0..pasos {
            f.paso(0.0, medio_periodo / pasos as Real);
        }
        assert!(
            (f.angulo_pendulo + 0.05).abs() < 1e-4,
            "angulo = {}",
            f.angulo_pendulo
        );
    }

    #[test]
    fn el_modelo_lineal_coincide_con_la_pendiente_del_no_lineal() {
        let mut f = furuta();
        f.friccion_brazo = 2e-3;
        f.friccion_pendulo = 1e-4;
        let (a, b) = f.lineal_arriba();
        let arriba = [0.0, PI, 0.0, 0.0];
        let paso = 1e-2;

        let parecidos = |lineal: Real, numerico: Real| {
            (lineal - numerico).abs() <= 0.02 * lineal.abs().max(1.0)
        };
        for columna in 0..4 {
            let (mut mas, mut menos) = (arriba, arriba);
            mas[columna] += paso;
            menos[columna] -= paso;
            let (mas, menos) = (f.derivada(&mas, 0.0), f.derivada(&menos, 0.0));
            for fila in 0..4 {
                let numerico = (mas[fila] - menos[fila]) / (2.0 * paso);
                let lineal = a.0[fila][columna];
                assert!(
                    parecidos(lineal, numerico),
                    "a[{fila}][{columna}]: {lineal} y {numerico}"
                );
            }
        }
        let (mas, menos) = (f.derivada(&arriba, paso), f.derivada(&arriba, -paso));
        for fila in 0..4 {
            let numerico = (mas[fila] - menos[fila]) / (2.0 * paso);
            assert!(parecidos(b.0[fila][0], numerico), "b[{fila}]: {numerico}");
        }
    }
}
