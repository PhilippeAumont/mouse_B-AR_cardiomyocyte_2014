# Model Information
This is the implementation in GNU Octave of the mouse cardiomyocyte model with Beta1-adrenoreceptor signaling from Bondarenko VE (2014). 
The model contains relevant channels/currents for sodium, potassium and calcium.Markov state models are used to model L-type Ca Channels (LCC), Fast Na channels, Ryanodine Receptors (RyR) and rapid delayed rectifier K channels (I Kr). 
This is also a compartmental model, using different compartments for the cytosol, diadic clefts, caveolae, and extracaveolae (which itself contains the diadic clefts). 
In addition, a model of lysosomal ion exchange is included. This adds two new compartments to the model: the lysosome and the microdomain, and small volume that represents the sum of all signalling microdomain found between lysosomes and the SR. 

## Use
The scripts `generate_parameters.m` and `set_starting_conditions.m` are used to generate files containing model parameters and initial state variable values respectively. The simulation settings can be set in the `main.m` file. The rest of the model is compartmentalized in different scripts in `+rates/` and `+odes/` and `+lys/`. Unless you want to change the model, you do not need to change anything is those directories. 

The model including the lysosomal ion transport is in the branch "Acidic-Ca-Storage".

## Model Sources
Bondarenko VE (2014) A Compartmentalized Mathematical Model of the β1-Adrenergic Signaling System in Mouse Ventricular Myocytes. PLoS ONE 9(2): e89113. https://doi.org/10.1371/journal.pone.0089113 

Astaburuaga, R., Quintanar Haro, O.D., Stauber, T., Relógio, A. A Mathematical Model of Lysosomal Ion Homeostasis Points to Differential Effects of Cl− Transport in Ca2+ Dynamics. Cells 2019, 8, 1263. https://doi.org/10.3390/cells8101263
