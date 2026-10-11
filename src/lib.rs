//! Bloques de control para sistemas embebidos.
//!
//! Todo funciona sin `std` y sin memoria dinámica: los tamaños se fijan al
//! compilar. Los cálculos van en [`Real`], que es `f32` salvo que se active la
//! característica `f64`.
//!
//! El tiempo se pasa siempre en segundos y cada bloque recibe su `dt`, así que
//! no dependen de ningún reloj ni de ninguna plataforma.

#![cfg_attr(not(test), no_std)]

pub mod estados;
pub mod filtro;
pub mod kalman;
pub mod matriz;
pub mod perfil;
pub mod pid;
pub mod planta;

#[cfg(not(feature = "f64"))]
pub type Real = f32;
#[cfg(feature = "f64")]
pub type Real = f64;

pub(crate) type Mate = libm::Libm<Real>;
pub(crate) const PI: Real = core::f64::consts::PI as Real;

// Compila los ejemplos del README junto con las pruebas.
#[cfg(doctest)]
#[doc = include_str!("../README.md")]
struct Leeme;
