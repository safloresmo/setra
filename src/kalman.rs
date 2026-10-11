//! Filtro de Kalman lineal.

use crate::matriz::{Matriz, Vector};

/// Estimador para el sistema discreto
///
/// ```text
/// x[k+1] = A·x[k] + B·u[k] + ruido de proceso (covarianza Q)
/// z[k]   = H·x[k] + ruido de medida (covarianza R)
/// ```
///
/// con `N` estados, `M` medidas y `U` entradas.
#[derive(Clone, Copy, Debug)]
pub struct Kalman<const N: usize, const M: usize, const U: usize> {
    pub a: Matriz<N, N>,
    pub b: Matriz<N, U>,
    pub h: Matriz<M, N>,
    pub q: Matriz<N, N>,
    pub r: Matriz<M, M>,
    /// Estado estimado.
    pub x: Vector<N>,
    /// Covarianza del error de la estimación.
    pub p: Matriz<N, N>,
}

impl<const N: usize, const M: usize, const U: usize> Kalman<N, M, U> {
    /// Avanza la estimación un periodo con la entrada aplicada.
    pub fn predecir(&mut self, u: &Vector<U>) {
        self.x = self.a * self.x + self.b * *u;
        self.p = self.a * self.p * self.a.transpuesta() + self.q;
    }

    /// Corrige la estimación con una medida. Devuelve `false`, sin tocar
    /// nada, si la covarianza de la innovación no se puede invertir.
    pub fn corregir(&mut self, z: &Vector<M>) -> bool {
        let ht = self.h.transpuesta();
        let Some(s_inversa) = (self.h * self.p * ht + self.r).inversa() else {
            return false;
        };
        let k = self.p * ht * s_inversa;

        self.x = self.x + k * (*z - self.h * self.x);
        // Forma de Joseph: mantiene P simétrica aunque se acumule redondeo.
        let ikh = Matriz::identidad() - k * self.h;
        self.p = ikh * self.p * ikh.transpuesta() + k * self.r * k.transpuesta();
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Real;

    // Ruido uniforme en [-1, 1] que se repite igual en cada ejecución.
    struct Ruido(u32);

    impl Ruido {
        fn siguiente(&mut self) -> Real {
            self.0 = self.0.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
            (self.0 >> 8) as Real / (1u32 << 23) as Real - 1.0
        }
    }

    #[test]
    fn estima_la_velocidad_midiendo_solo_la_posicion() {
        const DT: Real = 0.01;
        const VELOCIDAD: Real = 3.0;
        const AMPLITUD: Real = 0.2;

        let mut filtro = Kalman::<2, 1, 1> {
            a: Matriz([[1.0, DT], [0.0, 1.0]]),
            b: Matriz::ceros(),
            h: Matriz([[1.0, 0.0]]),
            q: Matriz::diagonal([1e-6, 1e-5]),
            // Varianza de un ruido uniforme de esa amplitud.
            r: Matriz([[AMPLITUD * AMPLITUD / 3.0]]),
            x: Vector::ceros(),
            p: Matriz::diagonal([1.0, 100.0]),
        };

        let mut ruido = Ruido(7);
        let mut posicion = 0.0;
        let mut peor_error: Real = 0.0;
        for n in 0..1000 {
            posicion += VELOCIDAD * DT;
            filtro.predecir(&Vector::ceros());
            assert!(filtro.corregir(&Vector::columna([posicion + AMPLITUD * ruido.siguiente()])));
            if n >= 500 {
                peor_error = peor_error.max((filtro.x.0[1][0] - VELOCIDAD).abs());
            }
        }

        assert!(peor_error < 0.15, "error de velocidad = {peor_error}");
        assert!((filtro.x.0[0][0] - posicion).abs() < 0.05);
        // La incertidumbre baja muy por debajo de la del sensor.
        assert!(filtro.p.0[0][0] < filtro.r.0[0][0] / 5.0);
    }

    #[test]
    fn la_entrada_mueve_la_prediccion() {
        let mut filtro = Kalman::<1, 1, 1> {
            a: Matriz([[1.0]]),
            b: Matriz([[0.5]]),
            h: Matriz([[1.0]]),
            q: Matriz([[0.0]]),
            r: Matriz([[1.0]]),
            x: Vector::columna([2.0]),
            p: Matriz([[1.0]]),
        };
        filtro.predecir(&Vector::columna([4.0]));
        assert_eq!(filtro.x.0[0][0], 4.0);
    }

    #[test]
    fn sin_ruido_declarado_ni_incertidumbre_no_se_puede_corregir() {
        let mut filtro = Kalman::<1, 1, 1> {
            a: Matriz([[1.0]]),
            b: Matriz::ceros(),
            h: Matriz([[1.0]]),
            q: Matriz::ceros(),
            r: Matriz::ceros(),
            x: Vector::columna([2.0]),
            p: Matriz::ceros(),
        };
        assert!(!filtro.corregir(&Vector::columna([5.0])));
        assert_eq!(filtro.x.0[0][0], 2.0);
    }
}
