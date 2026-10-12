// Banco de identificación de motores.
//
// Recibe por el puerto serie un ensayo, lo ejecuta registrando el encoder a
// ritmo fijo y devuelve los datos. Lo maneja `banco/identificar.py`.
//
// Orden:      ENSAYO <periodo en µs> <n> <ciclo 1> <ms 1> ... <ciclo n> <ms n>
//             Cada ciclo va en milésimas del voltaje de la fuente, de -1000 a 1000.
// Respuesta:  DATOS <muestras>, una línea "ciclo,cuentas" por muestra y FIN.
//             Si la orden no vale, una línea que empieza por ERROR.

#include <Arduino.h>
#include <ESP32Encoder.h>

#include "config.h"

const uint8_t CANAL_A = 0;
const uint8_t CANAL_B = 1;
const uint32_t CICLO_MAX = (1u << BITS_PWM) - 1;

struct Tramo {
  int16_t milesimas;
  uint32_t muestras;
};

static ESP32Encoder encoder;
static Tramo tramos[TRAMOS_MAX];
static int16_t *ciclos;
static int32_t *cuentas;

static void mover(int16_t milesimas) {
  uint32_t ciclo = (uint32_t)abs(milesimas) * CICLO_MAX / 1000;
#ifdef DRIVER_PWM_Y_SENTIDO
  digitalWrite(PIN_MOTOR_IN1, milesimas > 0);
  digitalWrite(PIN_MOTOR_IN2, milesimas < 0);
  ledcWrite(CANAL_A, ciclo);
#else
  ledcWrite(CANAL_A, milesimas > 0 ? ciclo : 0);
  ledcWrite(CANAL_B, milesimas < 0 ? ciclo : 0);
#endif
}

// Lee la orden y devuelve cuántas muestras pide, o 0 si no vale.
static size_t leerEnsayo(const String &orden, uint32_t &periodo, size_t &cantidad) {
  int leidos = 0;
  int n = 0;
  const char *texto = orden.c_str();
  if (sscanf(texto, "ENSAYO %u %d%n", &periodo, &n, &leidos) != 2) {
    Serial.println("ERROR no entiendo la orden");
    return 0;
  }
  if (periodo < 200 || n < 1 || n > (int)TRAMOS_MAX) {
    Serial.println("ERROR periodo o cantidad de tramos fuera de rango");
    return 0;
  }

  size_t total = 0;
  texto += leidos;
  for (int i = 0; i < n; i++) {
    int milesimas, duracion;
    if (sscanf(texto, "%d %d%n", &milesimas, &duracion, &leidos) != 2 || abs(milesimas) > 1000 ||
        duracion < 1) {
      Serial.println("ERROR tramo incompleto o fuera de rango");
      return 0;
    }
    texto += leidos;
    tramos[i] = {(int16_t)milesimas, (uint32_t)((uint64_t)duracion * 1000 / periodo)};
    total += tramos[i].muestras;
  }
  if (total == 0 || total > MUESTRAS_MAX) {
    Serial.println("ERROR el ensayo no cabe en la memoria");
    return 0;
  }
  cantidad = n;
  return total;
}

static void ensayar(uint32_t periodo, size_t cantidad, size_t total) {
  encoder.clearCount();
  size_t k = 0;
  unsigned long siguiente = micros();

  for (size_t i = 0; i < cantidad; i++) {
    for (uint32_t j = 0; j < tramos[i].muestras; j++) {
      while ((long)(micros() - siguiente) < 0) {
      }
      siguiente += periodo;
      // Primero se lee y después se aplica el voltaje, que queda puesto
      // hasta la muestra siguiente.
      cuentas[k] = (int32_t)encoder.getCount();
      ciclos[k] = tramos[i].milesimas;
      mover(tramos[i].milesimas);
      k++;
    }
  }
  mover(0);

  Serial.printf("DATOS %u\n", (unsigned)total);
  for (size_t i = 0; i < total; i++) {
    Serial.printf("%d,%ld\n", ciclos[i], (long)cuentas[i]);
  }
  Serial.println("FIN");
}

void setup() {
  Serial.begin(115200);
  Serial.setTimeout(50);

  ciclos = (int16_t *)malloc(MUESTRAS_MAX * sizeof(int16_t));
  cuentas = (int32_t *)malloc(MUESTRAS_MAX * sizeof(int32_t));

  encoder.attachFullQuad(PIN_ENCODER_A, PIN_ENCODER_B);

  pinMode(PIN_MOTOR_IN1, OUTPUT);
  pinMode(PIN_MOTOR_IN2, OUTPUT);
  ledcSetup(CANAL_A, FRECUENCIA_PWM, BITS_PWM);
#ifdef DRIVER_PWM_Y_SENTIDO
  ledcAttachPin(PIN_MOTOR_PWM, CANAL_A);
#else
  ledcSetup(CANAL_B, FRECUENCIA_PWM, BITS_PWM);
  ledcAttachPin(PIN_MOTOR_IN1, CANAL_A);
  ledcAttachPin(PIN_MOTOR_IN2, CANAL_B);
#endif
  mover(0);
}

void loop() {
  if (!Serial.available()) {
    return;
  }
  String orden = Serial.readStringUntil('\n');
  orden.trim();
  if (orden.isEmpty()) {
    return;
  }
  if (ciclos == nullptr || cuentas == nullptr) {
    Serial.println("ERROR sin memoria para el ensayo");
    return;
  }

  uint32_t periodo = 0;
  size_t cantidad = 0;
  size_t total = leerEnsayo(orden, periodo, cantidad);
  if (total > 0) {
    ensayar(periodo, cantidad, total);
  }
}
