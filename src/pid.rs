//! Controlador PID con salida limitada.

use crate::Real;

/// PID en forma paralela: `u = kp·e + ki·∫e + kd·de/dt`.
///
/// La derivada se toma de la medida y no del error, para que un escalón en la
/// referencia no dé un golpe en la salida; para seguir una referencia que se
/// mueve de forma suave conviene [`con_derivada_del_error`](Self::con_derivada_del_error).
/// Cuando la salida satura, la integral deja de crecer en ese sentido.
#[derive(Clone, Copy, Debug)]
pub struct Pid {
    pub kp: Real,
    pub ki: Real,
    pub kd: Real,
    salida_min: Real,
    salida_max: Real,
    tau_derivada: Real,
    derivada_del_error: bool,
    // Guarda ya el término completo (ki·∫e): cambiar `ki` en marcha no da saltos.
    integral: Real,
    derivada: Real,
    // Última muestra de la señal que se deriva: el error o la medida negada.
    previa: Option<Real>,
}

impl Pid {
    pub const fn nuevo(kp: Real, ki: Real, kd: Real) -> Self {
        Self {
            kp,
            ki,
            kd,
            salida_min: Real::NEG_INFINITY,
            salida_max: Real::INFINITY,
            tau_derivada: 0.0,
            derivada_del_error: false,
            integral: 0.0,
            derivada: 0.0,
            previa: None,
        }
    }

    /// Limita la salida al intervalo `[min, max]`.
    pub const fn con_limites(mut self, min: Real, max: Real) -> Self {
        self.salida_min = min;
        self.salida_max = max;
        self
    }

    /// Filtra la derivada con un paso bajo de primer orden de constante `tau`
    /// (en segundos). Con cero no se filtra.
    pub const fn con_filtro_derivada(mut self, tau: Real) -> Self {
        self.tau_derivada = tau;
        self
    }

    /// Deriva el error en vez de la medida. Sigue mejor una referencia en
    /// movimiento, a cambio de un golpe en la salida si la referencia salta.
    pub const fn con_derivada_del_error(mut self) -> Self {
        self.derivada_del_error = true;
        self
    }

    pub fn actualizar(&mut self, referencia: Real, medida: Real, dt: Real) -> Real {
        self.actualizar_con_adelanto(referencia, medida, 0.0, dt)
    }

    /// Igual que [`actualizar`](Self::actualizar), sumando a la salida un
    /// término de adelanto calculado fuera (por ejemplo `kv·velocidad`).
    pub fn actualizar_con_adelanto(
        &mut self,
        referencia: Real,
        medida: Real,
        adelanto: Real,
        dt: Real,
    ) -> Real {
        let error = referencia - medida;

        if dt > 0.0 {
            let derivable = if self.derivada_del_error {
                error
            } else {
                -medida
            };
            if let Some(previa) = self.previa {
                let cruda = (derivable - previa) / dt;
                let alfa = dt / (self.tau_derivada + dt);
                self.derivada += alfa * (cruda - self.derivada);
            }
            self.previa = Some(derivable);
        }

        let fijo = self.kp * error + self.kd * self.derivada + adelanto;
        let integral = self.integral + self.ki * error * dt.max(0.0);
        let sin_limite = fijo + integral;

        let empuja_la_saturacion = (sin_limite > self.salida_max && error > 0.0)
            || (sin_limite < self.salida_min && error < 0.0);
        if !empuja_la_saturacion {
            self.integral = integral;
        }

        (fijo + self.integral).clamp(self.salida_min, self.salida_max)
    }

    /// Término integral acumulado, en unidades de la salida.
    pub fn integral(&self) -> Real {
        self.integral
    }

    /// Borra la integral y la memoria de la derivada.
    pub fn reiniciar(&mut self) {
        self.integral = 0.0;
        self.derivada = 0.0;
        self.previa = None;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const DT: Real = 0.001;

    #[test]
    fn solo_proporcional() {
        let mut pid = Pid::nuevo(2.0, 0.0, 0.0);
        assert_eq!(pid.actualizar(1.0, 0.25, DT), 1.5);
    }

    #[test]
    fn la_salida_respeta_los_limites() {
        let mut pid = Pid::nuevo(100.0, 0.0, 0.0).con_limites(-12.0, 12.0);
        assert_eq!(pid.actualizar(1.0, 0.0, DT), 12.0);
        assert_eq!(pid.actualizar(-1.0, 0.0, DT), -12.0);
    }

    #[test]
    fn la_integral_no_crece_saturada() {
        let mut pid = Pid::nuevo(1.0, 50.0, 0.0).con_limites(-1.0, 1.0);
        for _ in 0..10_000 {
            pid.actualizar(10.0, 0.0, DT);
        }
        assert_eq!(pid.integral(), 0.0);

        // Al invertirse el error la salida deja la saturación enseguida.
        let salida = pid.actualizar(0.0, 0.5, DT);
        assert!(salida < 0.0, "salida = {salida}");
    }

    #[test]
    fn la_integral_elimina_el_error_permanente() {
        // Planta de primer orden con ganancia 0.5: sin integral quedaría error.
        let mut pid = Pid::nuevo(1.0, 20.0, 0.0);
        let mut y: Real = 0.0;
        for _ in 0..5_000 {
            let u = pid.actualizar(1.0, y, DT);
            y += DT * (0.5 * u - y) / 0.05;
        }
        assert!((y - 1.0).abs() < 1e-3, "y = {y}");
    }

    #[test]
    fn un_escalon_en_la_referencia_no_golpea_la_derivada() {
        let mut pid = Pid::nuevo(0.0, 0.0, 1.0);
        pid.actualizar(0.0, 0.0, DT);
        assert_eq!(pid.actualizar(100.0, 0.0, DT), 0.0);
    }

    #[test]
    fn la_derivada_del_error_sigue_una_rampa_sin_frenarla() {
        // Referencia y medida suben juntas: el error no cambia.
        let mut sobre_medida = Pid::nuevo(0.0, 0.0, 1.0);
        let mut sobre_error = Pid::nuevo(0.0, 0.0, 1.0).con_derivada_del_error();
        let (mut a, mut b) = (0.0, 0.0);
        for n in 0..10 {
            let x = n as Real * 0.5;
            a = sobre_medida.actualizar(x, x, 0.5);
            b = sobre_error.actualizar(x, x, 0.5);
        }
        assert_eq!(a, -1.0);
        assert_eq!(b, 0.0);
    }

    #[test]
    fn la_derivada_se_opone_al_movimiento() {
        let mut pid = Pid::nuevo(0.0, 0.0, 0.5);
        pid.actualizar(0.0, 0.0, DT);
        let salida = pid.actualizar(0.0, 0.01, DT);
        assert!((salida + 5.0).abs() < 1e-3, "salida = {salida}");
    }

    #[test]
    fn reiniciar_borra_la_memoria() {
        let mut pid = Pid::nuevo(0.0, 10.0, 1.0);
        for _ in 0..100 {
            pid.actualizar(1.0, 0.0, DT);
        }
        pid.reiniciar();
        assert_eq!(pid.integral(), 0.0);
        assert_eq!(pid.actualizar(0.0, 5.0, 0.0), 0.0);
    }
}
