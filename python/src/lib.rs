//! Módulo de Python: las mismas clases del núcleo, calculando en `f64`.

use nucleo::perfil::{Estado, Limites};
use nucleo::{filtro, perfil, pid, planta};
use pyo3::prelude::*;

/// PID en forma paralela con salida limitada y protección de la integral.
///
/// La derivada se toma de la medida salvo que se pida `derivada_del_error`,
/// que sigue mejor una referencia en movimiento.
#[pyclass(module = "setra")]
struct Pid(pid::Pid);

#[pymethods]
impl Pid {
    #[new]
    #[pyo3(signature = (kp, ki=0.0, kd=0.0, *, salida_min=f64::NEG_INFINITY,
                        salida_max=f64::INFINITY, tau_derivada=0.0, derivada_del_error=false))]
    fn nuevo(
        kp: f64,
        ki: f64,
        kd: f64,
        salida_min: f64,
        salida_max: f64,
        tau_derivada: f64,
        derivada_del_error: bool,
    ) -> Self {
        let pid = pid::Pid::nuevo(kp, ki, kd)
            .con_limites(salida_min, salida_max)
            .con_filtro_derivada(tau_derivada);
        Self(if derivada_del_error {
            pid.con_derivada_del_error()
        } else {
            pid
        })
    }

    /// Calcula la salida para un ciclo de `dt` segundos. `adelanto` se suma
    /// tal cual a la salida antes de limitarla.
    #[pyo3(signature = (referencia, medida, dt, adelanto=0.0))]
    fn actualizar(&mut self, referencia: f64, medida: f64, dt: f64, adelanto: f64) -> f64 {
        self.0
            .actualizar_con_adelanto(referencia, medida, adelanto, dt)
    }

    /// Borra la integral y la memoria de la derivada.
    fn reiniciar(&mut self) {
        self.0.reiniciar();
    }

    /// Término integral acumulado, en unidades de la salida.
    #[getter]
    fn integral(&self) -> f64 {
        self.0.integral()
    }

    #[getter]
    fn kp(&self) -> f64 {
        self.0.kp
    }

    #[setter(kp)]
    fn fijar_kp(&mut self, valor: f64) {
        self.0.kp = valor;
    }

    #[getter]
    fn ki(&self) -> f64 {
        self.0.ki
    }

    #[setter(ki)]
    fn fijar_ki(&mut self, valor: f64) {
        self.0.ki = valor;
    }

    #[getter]
    fn kd(&self) -> f64 {
        self.0.kd
    }

    #[setter(kd)]
    fn fijar_kd(&mut self, valor: f64) {
        self.0.kd = valor;
    }
}

/// Paso bajo de primer orden con constante de tiempo `tau`, en segundos.
#[pyclass(module = "setra")]
struct PasoBajo(filtro::PasoBajo);

#[pymethods]
impl PasoBajo {
    #[new]
    fn nuevo(tau: f64) -> Self {
        Self(filtro::PasoBajo::nuevo(tau))
    }

    /// Filtro con la frecuencia de corte dada en hercios.
    #[staticmethod]
    fn con_corte(frecuencia: f64) -> Self {
        Self(filtro::PasoBajo::con_corte(frecuencia))
    }

    fn actualizar(&mut self, entrada: f64, dt: f64) -> f64 {
        self.0.actualizar(entrada, dt)
    }

    fn reiniciar(&mut self) {
        self.0.reiniciar();
    }
}

/// Filtro de segundo orden para una frecuencia de muestreo fija. Se crea con
/// los coeficientes normalizados o con `paso_bajo` y `rechaza_banda`.
#[pyclass(module = "setra")]
struct Bicuadratico(filtro::Bicuadratico);

#[pymethods]
impl Bicuadratico {
    #[new]
    fn nuevo(b: [f64; 3], a: [f64; 2]) -> Self {
        Self(filtro::Bicuadratico::nuevo(b, a))
    }

    /// Paso bajo; con el `q` por omisión la respuesta es de Butterworth.
    #[staticmethod]
    #[pyo3(signature = (frecuencia, muestreo, q=std::f64::consts::FRAC_1_SQRT_2))]
    fn paso_bajo(frecuencia: f64, muestreo: f64, q: f64) -> Self {
        Self(filtro::Bicuadratico::paso_bajo(frecuencia, muestreo, q))
    }

    /// Rechaza una banda alrededor de `frecuencia`, más estrecha con más `q`.
    #[staticmethod]
    fn rechaza_banda(frecuencia: f64, muestreo: f64, q: f64) -> Self {
        Self(filtro::Bicuadratico::rechaza_banda(frecuencia, muestreo, q))
    }

    fn actualizar(&mut self, entrada: f64) -> f64 {
        self.0.actualizar(entrada)
    }

    fn reiniciar(&mut self) {
        self.0.reiniciar();
    }
}

/// Limita el cambio de una señal a `pendiente` unidades por segundo.
#[pyclass(module = "setra")]
struct LimitadorPendiente(filtro::LimitadorPendiente);

#[pymethods]
impl LimitadorPendiente {
    #[new]
    fn nuevo(pendiente: f64) -> Self {
        Self(filtro::LimitadorPendiente::nuevo(pendiente))
    }

    fn actualizar(&mut self, entrada: f64, dt: f64) -> f64 {
        self.0.actualizar(entrada, dt)
    }

    fn reiniciar(&mut self) {
        self.0.reiniciar();
    }
}

/// Movimiento de `desde` a `hasta` con aceleración constante, sin pasar de
/// la velocidad ni de la aceleración máximas.
#[pyclass(module = "setra")]
struct PerfilTrapezoidal(perfil::PerfilTrapezoidal);

#[pymethods]
impl PerfilTrapezoidal {
    #[new]
    #[pyo3(signature = (velocidad, aceleracion, desde, hasta, *,
                        velocidad_inicial=0.0, velocidad_final=0.0))]
    fn nuevo(
        velocidad: f64,
        aceleracion: f64,
        desde: f64,
        hasta: f64,
        velocidad_inicial: f64,
        velocidad_final: f64,
    ) -> Self {
        Self(perfil::PerfilTrapezoidal::planear(
            Limites {
                velocidad,
                aceleracion,
            },
            Estado {
                posicion: desde,
                velocidad: velocidad_inicial,
            },
            Estado {
                posicion: hasta,
                velocidad: velocidad_final,
            },
        ))
    }

    /// Posición y velocidad de referencia a los `t` segundos de empezar.
    fn estado(&self, t: f64) -> (f64, f64) {
        let estado = self.0.estado(t);
        (estado.posicion, estado.velocidad)
    }

    fn terminado(&self, t: f64) -> bool {
        self.0.terminado(t)
    }

    /// Segundos que tarda el movimiento completo.
    #[getter]
    fn duracion(&self) -> f64 {
        self.0.duracion()
    }
}

/// Motor de corriente continua de imanes permanentes, en unidades del SI.
#[pyclass(module = "setra")]
struct MotorDc(planta::MotorDc);

#[pymethods]
impl MotorDc {
    #[new]
    #[pyo3(signature = (resistencia, inductancia, kt, ke, inercia, friccion=0.0))]
    fn nuevo(
        resistencia: f64,
        inductancia: f64,
        kt: f64,
        ke: f64,
        inercia: f64,
        friccion: f64,
    ) -> Self {
        Self(planta::MotorDc {
            resistencia,
            inductancia,
            kt,
            ke,
            inercia,
            friccion,
            corriente: 0.0,
            velocidad: 0.0,
            posicion: 0.0,
        })
    }

    /// Avanza `dt` segundos, repartidos en `subpasos` pasos de integración.
    /// Cada paso debe durar menos que `inductancia / resistencia`.
    #[pyo3(signature = (voltaje, dt, *, par_carga=0.0, subpasos=1))]
    fn paso(&mut self, voltaje: f64, dt: f64, par_carga: f64, subpasos: u32) {
        let subpasos = subpasos.max(1);
        for _ in 0..subpasos {
            self.0.paso(voltaje, par_carga, dt / f64::from(subpasos));
        }
    }

    #[getter]
    fn corriente(&self) -> f64 {
        self.0.corriente
    }

    #[getter]
    fn velocidad(&self) -> f64 {
        self.0.velocidad
    }

    #[getter]
    fn posicion(&self) -> f64 {
        self.0.posicion
    }

    #[getter]
    fn resistencia(&self) -> f64 {
        self.0.resistencia
    }

    #[getter]
    fn inductancia(&self) -> f64 {
        self.0.inductancia
    }

    #[getter]
    fn kt(&self) -> f64 {
        self.0.kt
    }

    #[getter]
    fn ke(&self) -> f64 {
        self.0.ke
    }

    #[getter]
    fn inercia(&self) -> f64 {
        self.0.inercia
    }

    #[getter]
    fn friccion(&self) -> f64 {
        self.0.friccion
    }
}

/// Péndulo simple con par en el eje. El ángulo se mide desde la vertical
/// hacia abajo: `math.pi` es la posición invertida.
#[pyclass(module = "setra")]
struct Pendulo(planta::Pendulo);

#[pymethods]
impl Pendulo {
    #[new]
    #[pyo3(signature = (masa, longitud, *, friccion=0.0, gravedad=9.81, angulo=0.0, velocidad=0.0))]
    fn nuevo(
        masa: f64,
        longitud: f64,
        friccion: f64,
        gravedad: f64,
        angulo: f64,
        velocidad: f64,
    ) -> Self {
        Self(planta::Pendulo {
            masa,
            longitud,
            friccion,
            gravedad,
            angulo,
            velocidad,
        })
    }

    /// Avanza `dt` segundos con el par dado, en `subpasos` pasos de integración.
    #[pyo3(signature = (par, dt, *, subpasos=1))]
    fn paso(&mut self, par: f64, dt: f64, subpasos: u32) {
        let subpasos = subpasos.max(1);
        for _ in 0..subpasos {
            self.0.paso(par, dt / f64::from(subpasos));
        }
    }

    /// Energía mecánica, con el cero en el punto más bajo.
    #[getter]
    fn energia(&self) -> f64 {
        self.0.energia()
    }

    #[getter]
    fn angulo(&self) -> f64 {
        self.0.angulo
    }

    #[getter]
    fn velocidad(&self) -> f64 {
        self.0.velocidad
    }
}

#[pymodule(name = "setra")]
fn modulo(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Pid>()?;
    m.add_class::<PasoBajo>()?;
    m.add_class::<Bicuadratico>()?;
    m.add_class::<LimitadorPendiente>()?;
    m.add_class::<PerfilTrapezoidal>()?;
    m.add_class::<MotorDc>()?;
    m.add_class::<Pendulo>()?;
    Ok(())
}
