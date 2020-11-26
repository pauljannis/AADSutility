/*
ABSORBANCE SORTING
2D sorting gates:
-Takes absolute peaks --> low absorbent peaks are now sortable for negative selections
-Size selection --> Fusions or satellites can be excluded, reducing false positive rate
Author: Paul Zurek (pjz26@cam.ac.uk)
*/

//SETUP
//Sorting variables
float thresh = 8.0;     //set this to voltage ~0.5 below baseline!
float sortVH = 5.0;    //High voltage threshold
float sortVL = 0.0;    //Low voltage threshold
unsigned long peakWH = 2000UL;    //High size threshold (in us)
unsigned long peakWL = 500UL;     //Low size threshold (in us)
unsigned long minDist = 1000UL;   //min distance to last peak to exclude doublets

//Electrode timing
int electrode_delay = 250;   //Pulse delay in us
int electrode_pulse = 5;     //Pulse width in ms





//MAIN
void setup() {
  pinMode(13, OUTPUT);
  digitalWrite(13, LOW);
  REG_ADC_MR = (REG_ADC_MR & 0xFFF0FFFF) | 0x00020000;   //10x faster analogread (4us)
  Serial.begin(115200);                                  //due to faster ADC startup time
  analogReadResolution(12);
}

//Read and convert the detector voltage
float readout() { 
  int sensorValue = analogRead(A0); 
  float result = sensorValue / 4096.0 * 10;  //Scale to 0-10V
  return result;
}


//loop variables
float v = 10.0;
float vo = 10.0;
float peakV = 10.0;
unsigned long t1 = 0UL;         
unsigned long t2 = 0UL; 
unsigned long tlast = 0UL;


void loop() {
vo = v;
v = readout();     //get signal from detector and convert

if ((vo > thresh) and (v < thresh)) {   //left event border
  t1 = micros();
}
  
if ((vo < thresh) and (v > thresh)) {   //right event border
  tlast = t2;
  t2 = micros();

  //check size
  if (((t1-tlast) > minDist) && ((t2-t1) > peakWL) && ((t2-t1) < peakWH)) { 
   //check voltage
   if ((peakV > sortVL) && (peakV < sortVH)) {
     delayMicroseconds(electrode_delay);
     digitalWrite(13, HIGH);
     delay(electrode_pulse);
     digitalWrite(13, LOW);
    }
 }
}

//Below thresh control
if (v < thresh) {
  if (v < peakV) {
    peakV = v;
  }
}
//Above thresh control
if (v > thresh) {
  peakV = 10.0;
}

}

