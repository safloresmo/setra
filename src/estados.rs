//! Control por espacio de estados: discretización del modelo, realimentación
//! de estados óptima y observador.

use crate::matriz::{Matriz, Vector};
use crate::{Mate, Real};

const TERMINOS_SERIE: usize = 16;
const ITERACIONES_MAX: usize = 100_000;

/// Pasa el modelo continuo `dx/dt = A·x + B·u` al discreto
/// `x[k+1] = Ad·x[k] + Bd·u[k]`, con la entrada constante durante cada
/// periodo `dt`.
///
/// Usa la serie de la exponencial sin escalado previo, de modo que `dt` debe
/// ser pequeño frente a la dinámica del sistema, como lo es cualquier periodo
/// de control razonable.
pub fn discretizar<const N: usize, const U: usize>(
    a: &Matriz<N, N>,
    b: &Matriz<N, U>,
    dt: Real,
) -> (Matriz<N, N>, Matriz<N, U>) {
    // termino = (A·dt)^k / k!
    let mut termino = Matriz::<N, N>::identidad();
    let mut ad = termino;
    let mut integral = termino.escalada(dt);
    for k in 1..=TERMINOS_SERIE {
        termino = (termino * *a).escalada(dt / k as Real);
        ad = ad + termino;
        integral = integral + termino.escalada(dt / (k + 1) as Real);
    }
    (ad, integral * *b)
}

/// Ganancia `K` del regulador lineal cuadrático para el sistema discreto
/// `x[k+1] = A·x[k] + B·u[k]`: la ley `u = -K·x` minimiza la suma de
/// `xᵀ·Q·x + uᵀ·R·u`.
///
/// Itera la ecuación de Riccati hasta que deja de cambiar. Devuelve `None` si
/// no converge, que es lo que pasa cuando el sistema no se puede estabilizar
/// con esas entradas.
pub fn lqr<const N: usize, const U: usize>(
    a: &Matriz<N, N>,
    b: &Matriz<N, U>,
    q: &Matriz<N, N>,
    r: &Matriz<U, U>,
) -> Option<Matriz<U, N>> {
    let tolerancia = Mate::sqrt(Real::EPSILON);
    let (at, bt) = (a.transpuesta(), b.transpuesta());

    let mut p = *q;
    for _ in 0..ITERACIONES_MAX {
        let k = (*r + bt * p * *b).inversa()? * bt * p * *a;
        let siguiente = at * p * (*a - *b * k) + *q;
        // El redondeo la va dejando asimétrica; se promedia con su transpuesta.
        let siguiente = (siguiente + siguiente.transpuesta()).escalada(0.5);

        let (mut cambio, mut escala): (Real, Real) = (0.0, 1.0);
        for (nuevo, previo) in siguiente.0.iter().flatten().zip(p.0.iter().flatten()) {
            if !nuevo.is_finite() {
                return None;
            }
            cambio = cambio.max((nuevo - previo).abs());
            escala = escala.max(nuevo.abs());
        }
        p = siguiente;
        if cambio <= tolerancia * escala {
            return Some(k);
        }
    }
    None
}

/// Ganancia `L` para un [`Observador`] del sistema `x[k+1] = A·x[k]`,
/// `y[k] = C·x[k]`. `Q` y `R` pesan cuánto se desconfía del modelo y de las
/// medidas: con más `Q` el observador corrige más rápido y deja pasar más
/// ruido.
pub fn ganancia_observador<const N: usize, const M: usize>(
    a: &Matriz<N, N>,
    c: &Matriz<M, N>,
    q: &Matriz<N, N>,
    r: &Matriz<M, M>,
) -> Option<Matriz<N, M>> {
    // Diseñar el observador es el mismo problema que el regulador, transpuesto.
    lqr(&a.transpuesta(), &c.transpuesta(), q, r).map(|k| k.transpuesta())
}

/// Observador de Luenberger: reconstruye el estado completo a partir de las
/// entradas aplicadas y de las salidas medidas, con ganancia fija.
///
/// ```text
/// x[k+1] = A·x[k] + B·u[k] + L·(y[k] - C·x[k])
/// ```
#[derive(Clone, Copy, Debug)]
pub struct Observador<const N: usize, const M: usize, const U: usize> {
    pub a: Matriz<N, N>,
    pub b: Matriz<N, U>,
    pub c: Matriz<M, N>,
    pub l: Matriz<N, M>,
    /// Estado estimado.
    pub x: Vector<N>,
}

impl<const N: usize, const M: usize, const U: usize> Observador<N, M, U> {
    /// Avanza un periodo con la entrada aplicada y la medida de ese instante.
    pub fn actualizar(&mut self, u: &Vector<U>, y: &Vector<M>) {
        let innovacion = *y - self.c * self.x;
        self.x = self.a * self.x + self.b * *u + self.l * innovacion;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const DT: Real = 0.01;

    // Doble integrador: posición y velocidad, con la aceleración de entrada.
    fn doble_integrador() -> (Matriz<2, 2>, Matriz<2, 1>) {
        discretizar(
            &Matriz([[0.0, 1.0], [0.0, 0.0]]),
            &Matriz([[0.0], [1.0]]),
            DT,
        )
    }

    #[test]
    fn discretizar_el_doble_integrador() {
        let (ad, bd) = doble_integrador();
        assert_eq!(ad, Matriz([[1.0, DT], [0.0, 1.0]]));
        assert!((bd.0[0][0] - DT * DT / 2.0).abs() < 1e-9);
        assert!((bd.0[1][0] - DT).abs() < 1e-9);
    }

    #[test]
    fn discretizar_un_primer_orden() {
        let (ad, bd) = discretizar(&Matriz([[-2.0]]), &Matriz([[3.0]]), 0.1);
        let esperado = Mate::exp(-0.2);
        assert!((ad.0[0][0] - esperado).abs() < 1e-6);
        assert!((bd.0[0][0] - 3.0 * (1.0 - esperado) / 2.0).abs() < 1e-6);
    }

    #[test]
    fn lqr_de_un_integrador() {
        // Con todo en 1 la solución de Riccati es el número áureo.
        let uno = Matriz([[1.0]]);
        let k = lqr(&uno, &uno, &uno, &uno).unwrap();
        let aureo = (1.0 + Mate::sqrt(5.0)) / 2.0;
        assert!((k.0[0][0] - aureo / (1.0 + aureo)).abs() < 1e-3);
    }

    #[test]
    fn lqr_lleva_el_doble_integrador_al_origen() {
        let (a, b) = doble_integrador();
        let k = lqr(&a, &b, &Matriz::diagonal([10.0, 1.0]), &Matriz([[0.1]])).unwrap();

        let mut x = Vector::columna([1.0, 0.0]);
        for _ in 0..500 {
            x = (a - b * k) * x;
        }
        assert!(
            x.0[0][0].abs() < 1e-3 && x.0[1][0].abs() < 1e-3,
            "x = {x:?}"
        );
    }

    #[test]
    fn lqr_sin_solucion() {
        // Inestable y sin que la entrada llegue al estado.
        let k = lqr(
            &Matriz([[2.0]]),
            &Matriz([[0.0]]),
            &Matriz([[1.0]]),
            &Matriz([[1.0]]),
        );
        assert!(k.is_none());
    }

    #[test]
    fn el_observador_recupera_la_velocidad() {
        let (a, b) = doble_integrador();
        let c = Matriz([[1.0, 0.0]]);
        let l = ganancia_observador(&a, &c, &Matriz::diagonal([1.0, 100.0]), &Matriz([[0.01]]))
            .unwrap();
        let mut observador = Observador {
            a,
            b,
            c,
            l,
            x: Vector::ceros(),
        };

        // El sistema real arranca en otro estado y recibe una aceleración.
        let mut x = Vector::columna([2.0, -3.0]);
        let u = Vector::columna([0.5]);
        for _ in 0..300 {
            observador.actualizar(&u, &(c * x));
            x = a * x + b * u;
        }
        let error = x - observador.x;
        assert!(
            error.0[0][0].abs() < 1e-3 && error.0[1][0].abs() < 1e-3,
            "error = {error:?}"
        );
    }
}
