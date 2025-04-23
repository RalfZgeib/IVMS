;IVMS viscosity continuous measurements
;(Fixed extrusion times per feedrate)

;Use 3.1mL + 200s in MATLAB

;Measurement sequence (F>750-500-250-125-50-25)


T0

;Extrusion
G92 E0
G1 F750 E270
G4 S5
G1 F500 E450
G4 S5
G1 F250 E540
G4 S5
G1 F125 E585
G4 S5
G1 F50 E603
G4 S5
G1 F25 E612

;End

