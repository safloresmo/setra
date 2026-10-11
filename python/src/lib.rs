//! Módulo de Python: las mismas clases del núcleo, calculando en `f64`.

use nucleo::matriz::{Matriz, Vector};
use nucleo::perfil::{Estado, Limites, LimitesS};
use nucleo::{estados, filtro, kalman, perfil, pid, planta};
use pyo3::exceptions::PyValueError;
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

/// Movimiento de reposo a reposo de `desde` a `hasta` con el jerk limitado:
/// la aceleración sube y baja en rampa en vez de saltar.
#[pyclass(module = "setra")]
struct PerfilS(perfil::PerfilS);

#[pymethods]
impl PerfilS {
    #[new]
    fn nuevo(velocidad: f64, aceleracion: f64, jerk: f64, desde: f64, hasta: f64) -> Self {
        let limites = LimitesS {
            velocidad,
            aceleracion,
            jerk,
        };
        Self(perfil::PerfilS::planear(limites, desde, hasta))
    }

    /// Posición, velocidad y aceleración de referencia a los `t` segundos.
    fn estado(&self, t: f64) -> (f64, f64, f64) {
        let estado = self.0.estado(t);
        (estado.posicion, estado.velocidad, estado.aceleracion)
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

/// Péndulo invertido rotatorio, o de Furuta: un motor gira un brazo
/// horizontal que lleva en la punta un péndulo libre. El ángulo del péndulo
/// se mide desde la vertical hacia abajo: `math.pi` es la posición invertida.
#[pyclass(module = "setra")]
struct Furuta(planta::Furuta);

#[pymethods]
impl Furuta {
    #[new]
    #[pyo3(signature = (inercia_brazo, longitud_brazo, masa_pendulo, centro_pendulo, inercia_pendulo, *,
                        friccion_brazo=0.0, friccion_pendulo=0.0, gravedad=9.81,
                        angulo_brazo=0.0, angulo_pendulo=0.0))]
    #[allow(clippy::too_many_arguments)]
    fn nuevo(
        inercia_brazo: f64,
        longitud_brazo: f64,
        masa_pendulo: f64,
        centro_pendulo: f64,
        inercia_pendulo: f64,
        friccion_brazo: f64,
        friccion_pendulo: f64,
        gravedad: f64,
        angulo_brazo: f64,
        angulo_pendulo: f64,
    ) -> Self {
        Self(planta::Furuta {
            inercia_brazo,
            longitud_brazo,
            masa_pendulo,
            centro_pendulo,
            inercia_pendulo,
            friccion_brazo,
            friccion_pendulo,
            gravedad,
            angulo_brazo,
            angulo_pendulo,
            velocidad_brazo: 0.0,
            velocidad_pendulo: 0.0,
        })
    }

    /// Avanza `dt` segundos con el par dado en el eje del motor, en
    /// `subpasos` pasos de integración.
    #[pyo3(signature = (par, dt, *, subpasos=1))]
    fn paso(&mut self, par: f64, dt: f64, subpasos: u32) {
        let subpasos = subpasos.max(1);
        for _ in 0..subpasos {
            self.0.paso(par, dt / f64::from(subpasos));
        }
    }

    /// Modelo lineal continuo `(a, b)` alrededor del péndulo hacia arriba. El
    /// estado es el ángulo del brazo, la desviación del péndulo respecto a la
    /// vertical y sus dos velocidades; la entrada es el par.
    fn lineal_arriba(&self) -> (Filas, Filas) {
        let (a, b) = self.0.lineal_arriba();
        (filas(&a), filas(&b))
    }

    /// Energía mecánica, con el cero en el péndulo colgando y en reposo.
    #[getter]
    fn energia(&self) -> f64 {
        self.0.energia()
    }

    #[getter]
    fn angulo_brazo(&self) -> f64 {
        self.0.angulo_brazo
    }

    #[getter]
    fn angulo_pendulo(&self) -> f64 {
        self.0.angulo_pendulo
    }

    #[getter]
    fn velocidad_brazo(&self) -> f64 {
        self.0.velocidad_brazo
    }

    #[getter]
    fn velocidad_pendulo(&self) -> f64 {
        self.0.velocidad_pendulo
    }
}

// Los tamaños del núcleo se fijan al compilar. Una operación que depende de
// ellos implementa este rasgo y `segun_tamanos` la llama con los que tocan
// según las matrices recibidas.
trait PorTamanos {
    type Salida;

    fn llamar<const N: usize, const M: usize, const U: usize>(self) -> PyResult<Self::Salida>;
}

fn segun_tamanos<O: PorTamanos>(
    estados: usize,
    medidas: usize,
    entradas: usize,
    operacion: O,
) -> PyResult<O::Salida> {
    match estados {
        1 => segun_medidas::<O, 1>(medidas, entradas, operacion),
        2 => segun_medidas::<O, 2>(medidas, entradas, operacion),
        3 => segun_medidas::<O, 3>(medidas, entradas, operacion),
        4 => segun_medidas::<O, 4>(medidas, entradas, operacion),
        5 => segun_medidas::<O, 5>(medidas, entradas, operacion),
        6 => segun_medidas::<O, 6>(medidas, entradas, operacion),
        _ => Err(PyValueError::new_err("se admiten de 1 a 6 estados")),
    }
}

fn segun_medidas<O: PorTamanos, const N: usize>(
    medidas: usize,
    entradas: usize,
    operacion: O,
) -> PyResult<O::Salida> {
    match medidas {
        1 => segun_entradas::<O, N, 1>(entradas, operacion),
        2 => segun_entradas::<O, N, 2>(entradas, operacion),
        3 => segun_entradas::<O, N, 3>(entradas, operacion),
        4 => segun_entradas::<O, N, 4>(entradas, operacion),
        _ => Err(PyValueError::new_err("se admiten de 1 a 4 medidas")),
    }
}

fn segun_entradas<O: PorTamanos, const N: usize, const M: usize>(
    entradas: usize,
    operacion: O,
) -> PyResult<O::Salida> {
    match entradas {
        1 => operacion.llamar::<N, M, 1>(),
        2 => operacion.llamar::<N, M, 2>(),
        3 => operacion.llamar::<N, M, 3>(),
        4 => operacion.llamar::<N, M, 4>(),
        _ => Err(PyValueError::new_err("se admiten de 1 a 4 entradas")),
    }
}

type Filas = Vec<Vec<f64>>;

fn columnas(filas: &[Vec<f64>]) -> usize {
    filas.first().map_or(0, Vec::len)
}

fn matriz<const F: usize, const C: usize>(
    nombre: &str,
    filas: &[Vec<f64>],
) -> PyResult<Matriz<F, C>> {
    if filas.len() != F || filas.iter().any(|fila| fila.len() != C) {
        return Err(PyValueError::new_err(format!(
            "`{nombre}` debe ser de {F}×{C}"
        )));
    }
    let mut m = Matriz::ceros();
    for (destino, fila) in m.0.iter_mut().zip(filas) {
        destino.copy_from_slice(fila);
    }
    Ok(m)
}

fn vector<const N: usize>(nombre: &str, valores: &[f64]) -> PyResult<Vector<N>> {
    let valores: [f64; N] = valores
        .try_into()
        .map_err(|_| PyValueError::new_err(format!("`{nombre}` debe tener {N} valores")))?;
    Ok(Vector::columna(valores))
}

fn filas<const F: usize, const C: usize>(m: &Matriz<F, C>) -> Filas {
    m.0.iter().map(|fila| fila.to_vec()).collect()
}

fn valores<const N: usize>(v: &Vector<N>) -> Vec<f64> {
    v.0.iter().map(|fila| fila[0]).collect()
}

// Una entrada vacía vale por «sin entrada».
fn entrada<const U: usize>(u: &[f64]) -> PyResult<Vector<U>> {
    if u.is_empty() {
        Ok(Vector::ceros())
    } else {
        vector("u", u)
    }
}

struct Discretizar {
    a: Filas,
    b: Filas,
    dt: f64,
}

impl PorTamanos for Discretizar {
    type Salida = (Filas, Filas);

    fn llamar<const N: usize, const M: usize, const U: usize>(self) -> PyResult<Self::Salida> {
        let (ad, bd) =
            estados::discretizar::<N, U>(&matriz("a", &self.a)?, &matriz("b", &self.b)?, self.dt);
        Ok((filas(&ad), filas(&bd)))
    }
}

/// Pasa el modelo continuo `dx/dt = a·x + b·u` al discreto
/// `x[k+1] = ad·x[k] + bd·u[k]` con periodo `dt`. Devuelve `(ad, bd)`.
///
/// `dt` debe ser pequeño frente a la dinámica del sistema.
#[pyfunction]
fn discretizar(a: Filas, b: Filas, dt: f64) -> PyResult<(Filas, Filas)> {
    segun_tamanos(a.len(), 1, columnas(&b), Discretizar { a, b, dt })
}

struct Lqr {
    a: Filas,
    b: Filas,
    q: Filas,
    r: Filas,
}

impl PorTamanos for Lqr {
    type Salida = Filas;

    fn llamar<const N: usize, const M: usize, const U: usize>(self) -> PyResult<Self::Salida> {
        estados::lqr::<N, U>(
            &matriz("a", &self.a)?,
            &matriz("b", &self.b)?,
            &matriz("q", &self.q)?,
            &matriz("r", &self.r)?,
        )
        .map(|k| filas(&k))
        .ok_or_else(|| PyValueError::new_err("el regulador no converge con ese sistema"))
    }
}

/// Ganancia `k` del regulador lineal cuadrático para el sistema discreto
/// `x[k+1] = a·x[k] + b·u[k]`: la ley `u = -k·x` minimiza la suma de
/// `xᵀ·q·x + uᵀ·r·u`.
#[pyfunction]
fn lqr(a: Filas, b: Filas, q: Filas, r: Filas) -> PyResult<Filas> {
    segun_tamanos(a.len(), 1, columnas(&b), Lqr { a, b, q, r })
}

struct GananciaObservador {
    a: Filas,
    c: Filas,
    q: Filas,
    r: Filas,
}

impl PorTamanos for GananciaObservador {
    type Salida = Filas;

    fn llamar<const N: usize, const M: usize, const U: usize>(self) -> PyResult<Self::Salida> {
        estados::ganancia_observador::<N, M>(
            &matriz("a", &self.a)?,
            &matriz("c", &self.c)?,
            &matriz("q", &self.q)?,
            &matriz("r", &self.r)?,
        )
        .map(|l| filas(&l))
        .ok_or_else(|| PyValueError::new_err("el observador no converge con ese sistema"))
    }
}

/// Ganancia `l` para un `Observador` del sistema `x[k+1] = a·x[k]`,
/// `y[k] = c·x[k]`. Con más `q` frente a `r` el observador corrige más
/// rápido y deja pasar más ruido de la medida.
#[pyfunction]
fn ganancia_observador(a: Filas, c: Filas, q: Filas, r: Filas) -> PyResult<Filas> {
    segun_tamanos(a.len(), c.len(), 1, GananciaObservador { a, c, q, r })
}

trait Reconstructor: Send + Sync {
    fn actualizar(&mut self, u: &[f64], y: &[f64]) -> PyResult<()>;
    fn x(&self) -> Vec<f64>;
}

impl<const N: usize, const M: usize, const U: usize> Reconstructor
    for estados::Observador<N, M, U>
{
    fn actualizar(&mut self, u: &[f64], y: &[f64]) -> PyResult<()> {
        estados::Observador::actualizar(self, &entrada(u)?, &vector("y", y)?);
        Ok(())
    }

    fn x(&self) -> Vec<f64> {
        valores(&self.x)
    }
}

struct CrearObservador {
    a: Filas,
    b: Filas,
    c: Filas,
    l: Filas,
    x: Option<Vec<f64>>,
}

impl PorTamanos for CrearObservador {
    type Salida = Box<dyn Reconstructor>;

    fn llamar<const N: usize, const M: usize, const U: usize>(self) -> PyResult<Self::Salida> {
        Ok(Box::new(estados::Observador::<N, M, U> {
            a: matriz("a", &self.a)?,
            b: matriz("b", &self.b)?,
            c: matriz("c", &self.c)?,
            l: matriz("l", &self.l)?,
            x: match &self.x {
                Some(x) => vector("x", x)?,
                None => Vector::ceros(),
            },
        }))
    }
}

// El filtro de Kalman del núcleo, visto sin sus tamaños.
trait Estimador: Send + Sync {
    fn predecir(&mut self, u: &[f64]) -> PyResult<()>;
    fn corregir(&mut self, z: &[f64]) -> PyResult<bool>;
    fn x(&self) -> Vec<f64>;
    fn p(&self) -> Filas;
}

impl<const N: usize, const M: usize, const U: usize> Estimador for kalman::Kalman<N, M, U> {
    fn predecir(&mut self, u: &[f64]) -> PyResult<()> {
        kalman::Kalman::predecir(self, &entrada(u)?);
        Ok(())
    }

    fn corregir(&mut self, z: &[f64]) -> PyResult<bool> {
        Ok(kalman::Kalman::corregir(self, &vector("z", z)?))
    }

    fn x(&self) -> Vec<f64> {
        valores(&self.x)
    }

    fn p(&self) -> Filas {
        filas(&self.p)
    }
}

struct CrearKalman {
    a: Filas,
    b: Option<Filas>,
    h: Filas,
    q: Filas,
    r: Filas,
    x: Option<Vec<f64>>,
    p: Option<Filas>,
}

impl PorTamanos for CrearKalman {
    type Salida = Box<dyn Estimador>;

    fn llamar<const N: usize, const M: usize, const U: usize>(self) -> PyResult<Self::Salida> {
        Ok(Box::new(kalman::Kalman::<N, M, U> {
            a: matriz("a", &self.a)?,
            b: match &self.b {
                Some(b) => matriz("b", b)?,
                None => Matriz::ceros(),
            },
            h: matriz("h", &self.h)?,
            q: matriz("q", &self.q)?,
            r: matriz("r", &self.r)?,
            x: match &self.x {
                Some(x) => vector("x", x)?,
                None => Vector::ceros(),
            },
            p: match &self.p {
                Some(p) => matriz("p", p)?,
                None => Matriz::identidad(),
            },
        }))
    }
}

#[derive(FromPyObject)]
enum Valores {
    Uno(f64),
    Varios(Vec<f64>),
}

impl Valores {
    fn lista(self) -> Vec<f64> {
        match self {
            Self::Uno(valor) => vec![valor],
            Self::Varios(valores) => valores,
        }
    }
}

/// Filtro de Kalman lineal para el sistema discreto
///
///     x[k+1] = a·x[k] + b·u[k] + ruido de proceso (covarianza q)
///     z[k]   = h·x[k] + ruido de medida (covarianza r)
///
/// Las matrices se dan como listas de filas o arreglos de NumPy. Admite hasta
/// 6 estados, 4 medidas y 4 entradas. Sin `x` el estado inicial es cero y sin
/// `p` la covarianza inicial es la identidad.
#[pyclass(module = "setra")]
struct Kalman(Box<dyn Estimador>);

#[pymethods]
impl Kalman {
    #[new]
    #[pyo3(signature = (a, h, q, r, *, b=None, x=None, p=None))]
    fn nuevo(
        a: Filas,
        h: Filas,
        q: Filas,
        r: Filas,
        b: Option<Filas>,
        x: Option<Vec<f64>>,
        p: Option<Filas>,
    ) -> PyResult<Self> {
        // Sin `b` no hay entradas: queda una sola, con ganancia cero.
        let entradas = b.as_deref().map_or(1, columnas);
        let crear = CrearKalman {
            a,
            b,
            h,
            q,
            r,
            x,
            p,
        };
        segun_tamanos(crear.a.len(), crear.h.len(), entradas, crear).map(Self)
    }

    /// Avanza la estimación un periodo. Sin `u` la entrada es cero.
    #[pyo3(signature = (u=None))]
    fn predecir(&mut self, u: Option<Valores>) -> PyResult<()> {
        self.0.predecir(&u.map(Valores::lista).unwrap_or_default())
    }

    /// Corrige la estimación con una medida. Devuelve `False`, sin cambiar
    /// nada, si la covarianza de la innovación no se puede invertir.
    fn corregir(&mut self, z: Valores) -> PyResult<bool> {
        self.0.corregir(&z.lista())
    }

    /// Estado estimado.
    #[getter]
    fn x(&self) -> Vec<f64> {
        self.0.x()
    }

    /// Covarianza del error de la estimación, por filas.
    #[getter]
    fn p(&self) -> Filas {
        self.0.p()
    }
}

/// Observador de Luenberger: reconstruye el estado completo a partir de las
/// entradas aplicadas y de las salidas medidas, con ganancia fija `l`.
///
///     x[k+1] = a·x[k] + b·u[k] + l·(y[k] - c·x[k])
///
/// Admite los mismos tamaños que `Kalman`. Sin `x` el estado inicial es cero.
#[pyclass(module = "setra")]
struct Observador(Box<dyn Reconstructor>);

#[pymethods]
impl Observador {
    #[new]
    #[pyo3(signature = (a, b, c, l, *, x=None))]
    fn nuevo(a: Filas, b: Filas, c: Filas, l: Filas, x: Option<Vec<f64>>) -> PyResult<Self> {
        let crear = CrearObservador { a, b, c, l, x };
        segun_tamanos(crear.a.len(), crear.c.len(), columnas(&crear.b), crear).map(Self)
    }

    /// Avanza un periodo con la entrada aplicada y la medida de ese instante.
    fn actualizar(&mut self, u: Valores, y: Valores) -> PyResult<()> {
        self.0.actualizar(&u.lista(), &y.lista())
    }

    /// Estado estimado.
    #[getter]
    fn x(&self) -> Vec<f64> {
        self.0.x()
    }
}

#[pymodule(name = "setra")]
fn modulo(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Pid>()?;
    m.add_class::<PasoBajo>()?;
    m.add_class::<Bicuadratico>()?;
    m.add_class::<LimitadorPendiente>()?;
    m.add_class::<PerfilTrapezoidal>()?;
    m.add_class::<PerfilS>()?;
    m.add_class::<Kalman>()?;
    m.add_class::<Observador>()?;
    m.add_function(wrap_pyfunction!(discretizar, m)?)?;
    m.add_function(wrap_pyfunction!(lqr, m)?)?;
    m.add_function(wrap_pyfunction!(ganancia_observador, m)?)?;
    m.add_class::<MotorDc>()?;
    m.add_class::<Pendulo>()?;
    m.add_class::<Furuta>()?;
    Ok(())
}
