#include <Adafruit_Sensor.h>
#include <Adafruit_FXOS8700.h>
#include <Adafruit_FXAS21002C.h>
#include <Wire.h>
#include "SdFat.h"
#include "RingBuf.h"

#define PIN_SDA0 18  // Wire0 on Teensy 4.0
#define PIN_SCL0 19  // Wire0 on Teensy 4.0
#define PIN_SDA1 17  // Wire1 on Teensy 4.0
#define PIN_SCL1 16  // Wire1 on Teensy 4.0
#define PIN_RESET 15

/** Macro for mg per LSB at +/- 8g sensitivity (1 LSB = 0.000976mg) */
#define ACCEL_MG_LSB_8G (0.000976F)

Adafruit_FXOS8700 accelmag1 = Adafruit_FXOS8700(0x8700A, 0x8700B, 0x1F, Wire);
Adafruit_FXAS21002C gyro1 = Adafruit_FXAS21002C(0x0021002C, 0x21, Wire);
sensors_event_t aevent1, wevent1, mevent1;

Adafruit_FXOS8700 accelmag2 = Adafruit_FXOS8700(0x8700A, 0x8700B, 0x1F, Wire1); // shank
Adafruit_FXAS21002C gyro2 = Adafruit_FXAS21002C(0x0021002C, 0x21, Wire1);
sensors_event_t aevent2, wevent2, mevent2;

#define SD_CONFIG  SdioConfig(FIFO_SDIO)
#define RING_BUF_CAPACITY 400*512
#define LOG_FILE_SIZE 1000000000  // 150,000,000 bytes.

SdFs sd;
FsFile file;

// RingBuf for File type FsFile.
RingBuf<FsFile, RING_BUF_CAPACITY> rb;

char file_name[80];

struct packet {
  unsigned long int t;
  int16_t a1[3];
  int16_t w1[3];
  int16_t m1[3];
  int16_t a2[3];
  int16_t w2[3];
  int16_t m2[3];
} data;  // 40 bytes


void error(const char msg[]) {
  digitalWrite(LED_BUILTIN, HIGH);
  for (;;) {
    Serial.println(msg);
    delay(1000);
  }
}

void init_imu()
{
  // soft reset
  delay(50);
  digitalWrite(PIN_RESET, LOW);
  delay(50);
  digitalWrite(PIN_RESET, HIGH);
  delay(50);

  // IMU-1 =========================
    if (!accelmag1.begin(ACCEL_RANGE_8G) || !gyro1.begin(GYRO_RANGE_2000DPS)) {
      error("Ooops, no IMU-1 detected ... Check your wiring!");
    }
    Wire.setClock(400000);

  // IMU-1 =========================
  if (!accelmag2.begin(ACCEL_RANGE_8G) || !gyro2.begin(GYRO_RANGE_2000DPS)) {
    error("Ooops, no IMU-2 detected ... Check your wiring!");
  }
  Wire1.setClock(400000);

  delay(500);
}


void setup()
{
  pinMode(LED_BUILTIN, OUTPUT);
  digitalWrite(LED_BUILTIN, HIGH);

  // Open serial communications and wait for port to open:
  Serial.begin(2000000);
  delay(100);

  //      while(true) {
  //        Serial.println(sizeof(data));
  //        delay(1000);
  //      }
  init_imu();

  Serial.print("Initializing SD card...");
  if (!sd.begin(SD_CONFIG)) {
    error("initialization failed!");
    return;
  }

  // create new named file
  int file_name_count = 0;
  sprintf(file_name, "data%03d", file_name_count);
  while (sd.exists(file_name))
    sprintf(file_name, "data%03d", ++file_name_count);
  Serial.print("initialization succeeded! \nSaving to ");
  Serial.println(file_name);

  // open the file.
  if (!file.open(file_name, O_WRITE | O_CREAT))
    error("Cannot open file!");
  //  if (!file.preAllocate(LOG_FILE_SIZE))
  //     error("preAllocate failed\n");
  rb.begin(&file);

}

void loop()
{
  unsigned long tnow = micros();

  static unsigned long t0 = micros();
  const unsigned int dt0 = round(1e6 / 400);
  if (tnow >= t0) {

    gyro1.getEvent(&wevent1);
    accelmag1.getEvent(&aevent1, &mevent1);

    gyro2.getEvent(&wevent2);
    accelmag2.getEvent(&aevent2, &mevent2);

    // format packet ==========================================
    data.t = tnow;
    data.a1[0] = accelmag1.accel_raw.x;
    data.a1[1] = accelmag1.accel_raw.y;
    data.a1[2] = accelmag1.accel_raw.z;
    data.w1[0] = gyro1.raw.x;
    data.w1[1] = gyro1.raw.y;
    data.w1[2] = gyro1.raw.z;
    data.m1[0] = accelmag1.mag_raw.x;
    data.m1[1] = accelmag1.mag_raw.y;
    data.m1[2] = accelmag1.mag_raw.z;
    data.a2[0] = accelmag2.accel_raw.x;
    data.a2[1] = accelmag2.accel_raw.y;
    data.a2[2] = accelmag2.accel_raw.z;
    data.w2[0] = gyro2.raw.x;
    data.w2[1] = gyro2.raw.y;
    data.w2[2] = gyro2.raw.z;
    data.m2[0] = accelmag2.mag_raw.x;
    data.m2[1] = accelmag2.mag_raw.y;
    data.m2[2] = accelmag2.mag_raw.z;


    uint8_t *datablock = (uint8_t*) &data;
    rb.write("sop", 3);
    rb.write((const uint8_t *)&data, sizeof(data));

    // sync ring buffer =======================================
    // sync buffer, 512 block is optimal
    const int bsize = 512;
    static unsigned int data_count = 0;
    data_count += sizeof(data);
    if (data_count >= bsize) {
      data_count = 0;
      rb.sync();
    }

    // check if subject is static ==============================
    const int16_t a_tol = 0.5 / (ACCEL_MG_LSB_8G * 9.81);
    const int max_repeat_count = 2e6 / dt0;
    static int16_t a_last;
    static int repeat_count = 0;
    if (abs(accelmag2.accel_raw.x - a_last) < a_tol)
      repeat_count++;
    else {
      a_last = accelmag2.accel_raw.x;
      repeat_count = 0;
    }
    bool subject_is_still = repeat_count > max_repeat_count;

    // backup file ===============================================
    static unsigned long t_last_save = millis();
    const unsigned long period_save = 30e3;
    if ((millis() - t_last_save > period_save) || subject_is_still) {
      Serial.println("Saved!");
      file.flush();
      t_last_save = millis();
    }

    // check for sensor failure ==================================
    const int FAULT_LIMIT = 100;
    static int16_t w1_last, a1_last, m1_last, w2_last, a2_last, m2_last;
    static uint8_t w1_faultcount = 0, a1_faultcount = 0, m1_faultcount = 0,
                   w2_faultcount = 0, a2_faultcount = 0, m2_faultcount = 0;

    w1_faultcount = (gyro1.raw.x == w1_last)           ? (w1_faultcount + 1) : 0;
    a1_faultcount = (accelmag1.accel_raw.x == a1_last) ? (a1_faultcount + 1) : 0;
    m1_faultcount = (accelmag1.mag_raw.x == m1_last)   ? (m1_faultcount + 1) : 0;
    w2_faultcount = (gyro2.raw.x == w2_last)           ? (w2_faultcount + 1) : 0;
    a2_faultcount = (accelmag2.accel_raw.x == a2_last) ? (a2_faultcount + 1) : 0;
    m2_faultcount = (accelmag2.mag_raw.x == m2_last)   ? (m2_faultcount + 1) : 0;

    w1_last = gyro1.raw.x;
    a1_last = accelmag1.accel_raw.x;
    m1_last = accelmag1.mag_raw.x;
    w2_last = gyro2.raw.x;
    a2_last = accelmag2.accel_raw.x;
    m2_last = accelmag2.mag_raw.x;

    if (w1_faultcount > FAULT_LIMIT || a1_faultcount > FAULT_LIMIT || m1_faultcount > FAULT_LIMIT ||
        w2_faultcount > FAULT_LIMIT || a2_faultcount > FAULT_LIMIT || m2_faultcount > FAULT_LIMIT) {
      init_imu();
    }

    if (subject_is_still)
      digitalWrite(LED_BUILTIN, (millis() % 1000) < 500);
    else
      digitalWrite(LED_BUILTIN, (millis() % 200) < 100);

    tnow = micros();
    t0 = tnow + dt0 - ((tnow - t0) % dt0);
  }
}
