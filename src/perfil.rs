//! Perfiles de movimiento.

use crate::{Mate, Real};

/// Velocidad y aceleración máximas del movimiento, ambas positivas.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Limites {
    pub velocidad: Real,
    pub aceleracion: Real,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Estado {
    pub posicion: Real,
    pub velocidad: Real,
}

/// Movimiento entre dos estados con aceleración constante: acelera, mantiene
/// la velocidad máxima y frena. Si la distancia es corta el tramo de velocidad
/// constante desaparece y el perfil queda triangular.
///
/// El plan se calcula una vez y [`estado`](Self::estado) se consulta en cada
/// ciclo. Los estados inicial y final pueden llevar velocidad, de modo que se
/// puede replanear a mitad de un movimiento.
#[derive(Clone, Copy, Debug)]
pub struct PerfilTrapezoidal {
    limites: Limites,
    // Los cálculos se hacen como si el movimiento fuera siempre hacia
    // posiciones mayores; `sentido` deshace el cambio al entregar el estado.
    sentido: Real,
    inicial: Estado,
    meta: Estado,
    fin_aceleracion: Real,
    fin_crucero: Real,
    fin_frenado: Real,
}

impl PerfilTrapezoidal {
    pub fn planear(limites: Limites, inicial: Estado, meta: Estado) -> Self {
        let sentido = if inicial.posicion > meta.posicion {
            -1.0
        } else {
            1.0
        };
        let mut inicial = orientar(inicial, sentido);
        let meta = orientar(meta, sentido);
        inicial.velocidad = inicial.velocidad.min(limites.velocidad);

        let Limites {
            velocidad,
            aceleracion,
        } = limites;

        // Si se parte o se llega con velocidad, el perfil es un trapecio
        // completo al que le falta un trozo al principio o al final.
        let recorte_inicial = inicial.velocidad / aceleracion;
        let recorte_final = meta.velocidad / aceleracion;
        let distancia_completa = recorte_inicial * recorte_inicial * aceleracion / 2.0
            + (meta.posicion - inicial.posicion)
            + recorte_final * recorte_final * aceleracion / 2.0;

        let mut tiempo_aceleracion = velocidad / aceleracion;
        let mut distancia_crucero =
            distancia_completa - tiempo_aceleracion * tiempo_aceleracion * aceleracion;
        if distancia_crucero < 0.0 {
            tiempo_aceleracion = Mate::sqrt(distancia_completa / aceleracion);
            distancia_crucero = 0.0;
        }

        let fin_aceleracion = tiempo_aceleracion - recorte_inicial;
        let fin_crucero = fin_aceleracion + distancia_crucero / velocidad;
        let fin_frenado = fin_crucero + tiempo_aceleracion - recorte_final;

        Self {
            limites,
            sentido,
            inicial,
            meta,
            fin_aceleracion,
            fin_crucero,
            fin_frenado,
        }
    }

    /// Segundos que tarda el movimiento completo.
    pub fn duracion(&self) -> Real {
        self.fin_frenado
    }

    pub fn terminado(&self, t: Real) -> bool {
        t >= self.fin_frenado
    }

    /// Posición y velocidad de referencia a los `t` segundos de empezar.
    pub fn estado(&self, t: Real) -> Estado {
        let Limites {
            velocidad,
            aceleracion,
        } = self.limites;
        let t = t.max(0.0);

        let estado = if t < self.fin_aceleracion {
            Estado {
                posicion: self.inicial.posicion
                    + (self.inicial.velocidad + t * aceleracion / 2.0) * t,
                velocidad: self.inicial.velocidad + t * aceleracion,
            }
        } else if t < self.fin_crucero {
            Estado {
                posicion: self.inicial.posicion
                    + (self.inicial.velocidad + self.fin_aceleracion * aceleracion / 2.0)
                        * self.fin_aceleracion
                    + velocidad * (t - self.fin_aceleracion),
                velocidad,
            }
        } else if t < self.fin_frenado {
            let restante = self.fin_frenado - t;
            Estado {
                posicion: self.meta.posicion
                    - (self.meta.velocidad + restante * aceleracion / 2.0) * restante,
                velocidad: self.meta.velocidad + restante * aceleracion,
            }
        } else {
            self.meta
        };

        orientar(estado, self.sentido)
    }
}

fn orientar(estado: Estado, sentido: Real) -> Estado {
    Estado {
        posicion: estado.posicion * sentido,
        velocidad: estado.velocidad * sentido,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const LIMITES: Limites = Limites {
        velocidad: 2.0,
        aceleracion: 4.0,
    };
    const DT: Real = 1e-3;

    fn en(posicion: Real) -> Estado {
        Estado {
            posicion,
            velocidad: 0.0,
        }
    }

    // Recorre el perfil comprobando que respeta los límites y que la posición
    // es la integral de la velocidad.
    fn recorrer(perfil: &PerfilTrapezoidal) {
        let pasos = (perfil.duracion() / DT) as usize + 10;
        let mut previo = perfil.estado(0.0);
        for n in 1..=pasos {
            let actual = perfil.estado(n as Real * DT);
            assert!(actual.velocidad.abs() <= LIMITES.velocidad + 1e-4);
            let aceleracion = (actual.velocidad - previo.velocidad) / DT;
            assert!(
                aceleracion.abs() <= LIMITES.aceleracion * 1.01,
                "a = {aceleracion}"
            );
            let avance = (actual.velocidad + previo.velocidad) / 2.0 * DT;
            assert!((actual.posicion - previo.posicion - avance).abs() < 1e-4);
            previo = actual;
        }
    }

    #[test]
    fn trapecio_completo() {
        let perfil = PerfilTrapezoidal::planear(LIMITES, en(0.0), en(5.0));
        // 0.5 s para acelerar, 0.5 s para frenar y 4 m a 2 m/s.
        assert!((perfil.duracion() - 3.0).abs() < 1e-5);
        assert_eq!(perfil.estado(1.5).velocidad, 2.0);
        assert_eq!(perfil.estado(10.0), en(5.0));
        recorrer(&perfil);
    }

    #[test]
    fn distancia_corta_da_un_triangulo() {
        let perfil = PerfilTrapezoidal::planear(LIMITES, en(0.0), en(0.25));
        assert!((perfil.duracion() - 0.5).abs() < 1e-5);
        // Solo alcanza 1 m/s, a mitad del recorrido.
        assert!((perfil.estado(0.25).velocidad - 1.0).abs() < 1e-5);
        recorrer(&perfil);
    }

    #[test]
    fn hacia_atras() {
        let perfil = PerfilTrapezoidal::planear(LIMITES, en(3.0), en(-2.0));
        assert!((perfil.duracion() - 3.0).abs() < 1e-5);
        assert_eq!(perfil.estado(1.5).velocidad, -2.0);
        assert_eq!(perfil.estado(3.0), en(-2.0));
        recorrer(&perfil);
    }

    #[test]
    fn replanear_en_marcha_no_da_saltos() {
        let primero = PerfilTrapezoidal::planear(LIMITES, en(0.0), en(5.0));
        let a_medias = primero.estado(1.0);
        let segundo = PerfilTrapezoidal::planear(LIMITES, a_medias, en(8.0));
        assert_eq!(segundo.estado(0.0), a_medias);
        assert_eq!(segundo.estado(segundo.duracion()), en(8.0));
        recorrer(&segundo);
    }

    #[test]
    fn sin_distancia_no_hay_movimiento() {
        let perfil = PerfilTrapezoidal::planear(LIMITES, en(1.0), en(1.0));
        assert_eq!(perfil.duracion(), 0.0);
        assert_eq!(perfil.estado(0.0), en(1.0));
        assert!(perfil.terminado(0.0));
    }
}
