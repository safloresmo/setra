//! Matrices de tamaño fijo, lo justo para los estimadores.

use core::ops::{Add, Mul, Sub};

use crate::Real;

/// Matriz de `F` filas y `C` columnas, guardada por filas.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Matriz<const F: usize, const C: usize>(pub [[Real; C]; F]);

pub type Vector<const N: usize> = Matriz<N, 1>;

impl<const F: usize, const C: usize> Matriz<F, C> {
    pub const fn ceros() -> Self {
        Self([[0.0; C]; F])
    }

    pub fn transpuesta(&self) -> Matriz<C, F> {
        let mut t = Matriz::ceros();
        for f in 0..F {
            for c in 0..C {
                t.0[c][f] = self.0[f][c];
            }
        }
        t
    }
}

impl<const N: usize> Matriz<N, N> {
    pub fn identidad() -> Self {
        Self::diagonal([1.0; N])
    }

    pub fn diagonal(valores: [Real; N]) -> Self {
        let mut m = Self::ceros();
        for (i, valor) in valores.into_iter().enumerate() {
            m.0[i][i] = valor;
        }
        m
    }

    /// Inversa por Gauss-Jordan con pivoteo parcial. `None` si es singular.
    pub fn inversa(&self) -> Option<Self> {
        let mut a = self.0;
        let mut inversa = Self::identidad().0;

        for col in 0..N {
            let mut pivote = col;
            for f in col + 1..N {
                if a[f][col].abs() > a[pivote][col].abs() {
                    pivote = f;
                }
            }
            let divisor = a[pivote][col];
            if divisor == 0.0 || !divisor.is_finite() {
                return None;
            }
            a.swap(col, pivote);
            inversa.swap(col, pivote);

            for j in 0..N {
                a[col][j] /= divisor;
                inversa[col][j] /= divisor;
            }
            for f in 0..N {
                let factor = a[f][col];
                if f == col || factor == 0.0 {
                    continue;
                }
                for j in 0..N {
                    a[f][j] -= factor * a[col][j];
                    inversa[f][j] -= factor * inversa[col][j];
                }
            }
        }
        Some(Self(inversa))
    }
}

impl<const N: usize> Vector<N> {
    pub fn columna(valores: [Real; N]) -> Self {
        let mut v = Self::ceros();
        for (fila, valor) in v.0.iter_mut().zip(valores) {
            fila[0] = valor;
        }
        v
    }
}

impl<const F: usize, const K: usize, const C: usize> Mul<Matriz<K, C>> for Matriz<F, K> {
    type Output = Matriz<F, C>;

    fn mul(self, otra: Matriz<K, C>) -> Matriz<F, C> {
        let mut m = Matriz::ceros();
        for f in 0..F {
            for c in 0..C {
                for k in 0..K {
                    m.0[f][c] += self.0[f][k] * otra.0[k][c];
                }
            }
        }
        m
    }
}

impl<const F: usize, const C: usize> Add for Matriz<F, C> {
    type Output = Self;

    fn add(mut self, otra: Self) -> Self {
        for f in 0..F {
            for c in 0..C {
                self.0[f][c] += otra.0[f][c];
            }
        }
        self
    }
}

impl<const F: usize, const C: usize> Sub for Matriz<F, C> {
    type Output = Self;

    fn sub(mut self, otra: Self) -> Self {
        for f in 0..F {
            for c in 0..C {
                self.0[f][c] -= otra.0[f][c];
            }
        }
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parecidas<const N: usize>(a: &Matriz<N, N>, b: &Matriz<N, N>) -> bool {
        (0..N).all(|f| (0..N).all(|c| (a.0[f][c] - b.0[f][c]).abs() < 1e-4))
    }

    #[test]
    fn producto_de_tamanos_distintos() {
        let a = Matriz([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]);
        let b = Matriz([[7.0, 8.0], [9.0, 10.0], [11.0, 12.0]]);
        assert_eq!(a * b, Matriz([[58.0, 64.0], [139.0, 154.0]]));
    }

    #[test]
    fn transpuesta() {
        let a = Matriz([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]);
        assert_eq!(
            a.transpuesta(),
            Matriz([[1.0, 4.0], [2.0, 5.0], [3.0, 6.0]])
        );
    }

    #[test]
    fn inversa_por_la_matriz_da_la_identidad() {
        // El cero en la esquina obliga a pivotear.
        let a = Matriz([[0.0, 2.0, 1.0], [1.0, -1.0, 3.0], [4.0, 0.5, -2.0]]);
        let inversa = a.inversa().unwrap();
        assert!(parecidas(&(a * inversa), &Matriz::identidad()));
    }

    #[test]
    fn una_matriz_singular_no_tiene_inversa() {
        let a = Matriz([[1.0, 2.0], [2.0, 4.0]]);
        assert!(a.inversa().is_none());
    }
}
