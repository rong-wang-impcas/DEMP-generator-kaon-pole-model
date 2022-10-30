# DEMP-generator-kaon-pole-model
Event generator for deep exclusive kaon production based on kaon pole formula and kaon form factor

Event generation of  e p --> e' Lambda K^{+}  data

# Usage
Using ROOT to execute the codes:
>root -l test.cpp

Using ACLiC to compile and execute
>root -l

[0] .x test.cpp+



# In this version (v1.0), the beam-crossing angle is implemented.

The electron goes to the z direction, consistent with the fixed-target experiment.
see PionExculsiveElectroproduction() or SetBeamCrossAngle(double _angle)
for details of how the crossing angle is implemented.



Download ROOT and get started: https://root.cern.ch/

coding issue contact: rwang@impcas.ac.cn