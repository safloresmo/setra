#pragma once

#include <Arduino.h>

// Encoder en cuadratura del motor.
const uint8_t PIN_ENCODER_A = 34;
const uint8_t PIN_ENCODER_B = 35;

// Driver del motor. Hay dos formas de mandarlo:
//
//  - Con una señal de PWM y dos de sentido (L298N, TB6612): deja definido
//    DRIVER_PWM_Y_SENTIDO.
//  - Con el PWM directamente en sus dos entradas (DRV8833, MX1508): comenta
//    esa línea; PIN_MOTOR_PWM no se usa.
#define DRIVER_PWM_Y_SENTIDO

const uint8_t PIN_MOTOR_PWM = 25;
const uint8_t PIN_MOTOR_IN1 = 26;
const uint8_t PIN_MOTOR_IN2 = 27;

const uint32_t FRECUENCIA_PWM = 20000; // Hz, por encima de lo audible
const uint8_t BITS_PWM = 10;

// Lo más largo que puede durar un ensayo, en muestras.
const size_t MUESTRAS_MAX = 12000;
const size_t TRAMOS_MAX = 32;
