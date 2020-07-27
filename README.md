# O2SpaceChargeCalc
repository for implementation of calculation of the space charge distortions and corrections. 
The current implementation of the global distortions and corrections are performed in the class AliTPCSpaceCharge3DCalc.cxx in O2. However, this class is poorly optimized and contains several inconsistencies between the calculation of the distortions and corrections. The aim of this repository is to create a newly written, fast and precise calculation of the global distortions and corrections.
