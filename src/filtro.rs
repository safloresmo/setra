//! Filtros para señales muestreadas.

use crate::{Mate, PI, Real};

/// Paso bajo de primer orden. Admite un `dt` distinto en cada muestra.
#[derive(Clone, Copy, Debug)]
pub struct PasoBajo {
    tau: Real,
    salida: Option<Real>,
}

impl PasoBajo {
    /// `tau` es la constante de tiempo en segundos.
    pub const fn nuevo(tau: Real) -> Self {
        Self { tau, salida: None }
    }

    /// Filtro con la frecuencia de corte dada en hercios.
    pub fn con_corte(frecuencia: Real) -> Self {
        Self::nuevo(1.0 / (2.0 * PI * frecuencia))
    }

    /// La primera muestra pasa tal cual, para no arrancar desde cero.
    pub fn actualizar(&mut self, entrada: Real, dt: Real) -> Real {
        let salida = match self.salida {
            Some(previa) => previa + dt / (self.tau + dt) * (entrada - previa),
            None => entrada,
        };
        self.salida = Some(salida);
        salida
    }

    pub fn reiniciar(&mut self) {
        self.salida = None;
    }
}

/// Filtro de segundo orden para una frecuencia de muestreo fija.
#[derive(Clone, Copy, Debug)]
pub struct Bicuadratico {
    b: [Real; 3],
    a: [Real; 2],
    estado: [Real; 2],
}

impl Bicuadratico {
    /// Coeficientes ya normalizados (`a0 = 1`).
    pub const fn nuevo(b: [Real; 3], a: [Real; 2]) -> Self {
        Self {
            b,
            a,
            estado: [0.0; 2],
        }
    }

    /// Paso bajo con corte en `frecuencia`. Con `q = 0.7071` la respuesta es
    /// de Butterworth.
    pub fn paso_bajo(frecuencia: Real, muestreo: Real, q: Real) -> Self {
        let (coseno, alfa) = Self::parametros(frecuencia, muestreo, q);
        let a0 = 1.0 + alfa;
        let b1 = (1.0 - coseno) / a0;
        Self::nuevo(
            [b1 / 2.0, b1, b1 / 2.0],
            [-2.0 * coseno / a0, (1.0 - alfa) / a0],
        )
    }

    /// Rechaza una banda estrecha alrededor de `frecuencia`, más estrecha
    /// cuanto mayor es `q`. Sirve para quitar una resonancia mecánica.
    pub fn rechaza_banda(frecuencia: Real, muestreo: Real, q: Real) -> Self {
        let (coseno, alfa) = Self::parametros(frecuencia, muestreo, q);
        let a0 = 1.0 + alfa;
        let b1 = -2.0 * coseno / a0;
        Self::nuevo([1.0 / a0, b1, 1.0 / a0], [b1, (1.0 - alfa) / a0])
    }

    fn parametros(frecuencia: Real, muestreo: Real, q: Real) -> (Real, Real) {
        let w0 = 2.0 * PI * frecuencia / muestreo;
        (Mate::cos(w0), Mate::sin(w0) / (2.0 * q))
    }

    pub fn actualizar(&mut self, entrada: Real) -> Real {
        let salida = self.b[0] * entrada + self.estado[0];
        self.estado[0] = self.b[1] * entrada - self.a[0] * salida + self.estado[1];
        self.estado[1] = self.b[2] * entrada - self.a[1] * salida;
        salida
    }

    pub fn reiniciar(&mut self) {
        self.estado = [0.0; 2];
    }
}

/// Limita la rapidez con la que cambia una señal.
#[derive(Clone, Copy, Debug)]
pub struct LimitadorPendiente {
    pendiente: Real,
    salida: Option<Real>,
}

impl LimitadorPendiente {
    /// `pendiente` es el cambio máximo por segundo, en cualquier sentido.
    pub const fn nuevo(pendiente: Real) -> Self {
        Self {
            pendiente,
            salida: None,
        }
    }

    pub fn actualizar(&mut self, entrada: Real, dt: Real) -> Real {
        let salida = match self.salida {
            Some(previa) => {
                let paso = self.pendiente * dt;
                previa + (entrada - previa).clamp(-paso, paso)
            }
            None => entrada,
        };
        self.salida = Some(salida);
        salida
    }

    pub fn reiniciar(&mut self) {
        self.salida = None;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const MUESTREO: Real = 1000.0;
    const BUTTERWORTH: Real = core::f64::consts::FRAC_1_SQRT_2 as Real;

    // Amplitud de la salida ante un seno, medida ya pasado el transitorio.
    // Sale del valor eficaz: el máximo de las muestras se queda corto cuando
    // ninguna cae justo en la cresta.
    fn amplitud(filtro: &mut Bicuadratico, frecuencia: Real) -> Real {
        let mut suma = 0.0;
        for n in 0..4000 {
            let t = n as Real / MUESTREO;
            let salida = filtro.actualizar(Mate::sin(2.0 * PI * frecuencia * t));
            if n >= 2000 {
                suma += salida * salida;
            }
        }
        Mate::sqrt(2.0 * suma / 2000.0)
    }

    #[test]
    fn paso_bajo_alcanza_el_63_por_ciento_en_tau() {
        let mut filtro = PasoBajo::nuevo(0.1);
        filtro.actualizar(0.0, 1e-4);
        let mut salida = 0.0;
        for _ in 0..1000 {
            salida = filtro.actualizar(1.0, 1e-4);
        }
        assert!((salida - 0.632).abs() < 0.01, "salida = {salida}");
    }

    #[test]
    fn paso_bajo_arranca_en_la_primera_muestra() {
        let mut filtro = PasoBajo::con_corte(10.0);
        assert_eq!(filtro.actualizar(3.0, 0.001), 3.0);
    }

    #[test]
    fn bicuadratico_paso_bajo() {
        let nuevo = || Bicuadratico::paso_bajo(50.0, MUESTREO, BUTTERWORTH);
        assert!((amplitud(&mut nuevo(), 5.0) - 1.0).abs() < 0.01);
        // En el corte un Butterworth cae 3 dB.
        assert!((amplitud(&mut nuevo(), 50.0) - BUTTERWORTH).abs() < 0.02);
        assert!(amplitud(&mut nuevo(), 400.0) < 0.02);
    }

    #[test]
    fn bicuadratico_rechaza_banda() {
        let nuevo = || Bicuadratico::rechaza_banda(60.0, MUESTREO, 5.0);
        assert!(amplitud(&mut nuevo(), 60.0) < 0.01);
        assert!((amplitud(&mut nuevo(), 5.0) - 1.0).abs() < 0.02);
        assert!((amplitud(&mut nuevo(), 300.0) - 1.0).abs() < 0.02);
    }

    #[test]
    fn limitador_de_pendiente() {
        let mut limitador = LimitadorPendiente::nuevo(2.0);
        assert_eq!(limitador.actualizar(0.0, 0.1), 0.0);
        assert!((limitador.actualizar(10.0, 0.1) - 0.2).abs() < 1e-6);
        assert!((limitador.actualizar(10.0, 0.1) - 0.4).abs() < 1e-6);
        assert!((limitador.actualizar(0.3, 0.1) - 0.3).abs() < 1e-6);
    }
}
