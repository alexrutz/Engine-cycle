# engine-cycle
This is a student project to calculate the thermodynamic cycle of a an aircraft engine. Additionally, preliminary designs for compressor/turbine blades can be obtained. Work in Progress.

<p align="center">
  <img width="450" height="450" src="https://github.com/alexrutz/engine-cycle/blob/main/images/output_cycle_diagramm_T_s.png">
</p>

All you need to get this project running is Python, PyQt5, SciPy and matplotlib. First, download and install Python from https://www.python.org/downloads/. Open a terminal and open the directory where the engine-cycle source code lies. Type:
```
pip install -r requirements.txt
```
This will install all the needed packages and depencies in one go.\
If you want to manually install the packages, you can just install the three packages mentioned above with
```
pip install -package name-
```
This will also provide the needed dependencies.\
Now all you need to do is double click "main.py" and the main window should open. Alternatively, you can open "main.py" inside your IDE and hit run.

# Disclaimer
The preliminary design area for compressor and turbine blades is currently not finished yet and can produce false results. The main table for the thermodynamic cycle and the diagrams are accurate. The limits for input data as to where feasable results are produced have not been evaluated, but if you use common sense and dont accelerate the aircraft to Mach 12 it should be fine ;)

# Structure
This programm has a simple structure. The two main components are the input area and the result area. Between these two lies the main code. The "Magie"-button is used to process all the input data with formulas and give out the results.

<p align="center">
  <img width="740" height="450" src="https://github.com/alexrutz/engine-cycle/blob/main/images/overview.png">
</p>

# Input Area

This is where you put all the data that you have. If you are missing data for a specific input, you can just leave the default value (stems from GE CF34-81), which will probably be a reasonable gap filler.

### thermodynamic input

<p align="center">
  <img width="740" height="146" src="https://github.com/alexrutz/engine-cycle/blob/main/images/input_thermo.png">
</p>

"Standardtag" just defines the ambient conditions. Since most of the time you define altitude and Mach-number, this option is not being processed in the current state. Can be adjusted in the source code.\
Tamb: ambient temperature\
pamp: ambient pressure

"Betriebspunkt" operating point.\
M0: Mach-number of the aircraft\
H0: Altitude of the aircraft

"Druckaufteilung Kern" sets the pressure relation between low- and high-pressure-compressor as well as the relation between stage one and stage two of the high-pressure-compressor.\
VN / VH:   pressure relation between low- and high-pressure-compressor\
VH1 / VH2: pressure relation between state one and stage two of the high-pressure-compressor

"Nebenstromdaten" describes the fan operation\
lambda: bypass-ratio\
PI_of:  pressure ratio of the fan

"Luft/Kerosen": fuel and air characteristics\
Tref:	reference Temperature for air properties\
pref:	reference pressure for air properties\
RG:		specific gas constant air\
HU:		calorific value of fuel\
OTDF:	outlet-temperature-distribution-factor

"Prozess" overall cycle parameters\
Tt4:   burner exit temperature\
PI_K:  overall pressure ration\
Schub: thrust\
mp2:   inlet mass flow

"Eintrittsmassestrom": this decides whether mp2 is considered core mass flow or the entire massflow that goes through the fan\
Kernstrom:   core mass flow\
gesamtes TW: through fan

"Wirkungsgrade polytrop" polytropic efficiency of components\
F:     fan\
VN:    low-pressure-compressor\
VH     high-pressure-compressor (identical for both stages)\
TH:    high-pressure-turbine\
TN:    low-pressure-turbine\
BK:    burner\
Welle: mechanical shaft efficiency (bearings loss, ...)

"Druckverluste" total pressure loss over stations in engine\
BK: burner\
PDK: nozzle core\
PDN: nozzle bypass-ratio\
Einlauf: inlet\
mV1: cooling air mixture behind burner\
mHPT: cooling air mixture behind high-pressure-turbine\
mLPT: cooling air mixture behind low-pressure-turbine

"Kühlluft" cooling air ratio from mp2\
y_401: cooling air mixture behind burner\
y_41: cooling air mixture behind high-pressure-turbine\
y_5: cooling air mixture behind low-pressure-turbine

"Düse" this decides whether the nozzle is adapted (expanding to ambient pressure) or not\
angepasst:       adapted\
nicht angepasst: not adapted