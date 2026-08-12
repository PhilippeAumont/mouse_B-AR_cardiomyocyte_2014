# mouse_B-AR_cardiomyocyte
This is the implementation in GNU Octave of the mouse cardiomyocyte model with Beta1-adrenoreceptor signaling from Bondarenko VE (2014). 
The model contains relevant channels/currents for sodium, potassium and calcium.Markov state models are used to model L-type Ca Channels (LCC), Fast Na channels, Ryanodine Receptors (RyR) and rapid delayed rectifier K channels (I Kr). 
This is also a compartmental model, using different compartments for the cytosol, diadic clefts, caveolae, and extracaveolae (which itself contains the diadic clefts). 

## Use
The main.m script contains a set of parameters used to set up the simulations. The rest of the model is compartmentalized in different scripts in `+rates/` and `+odes/`.

## Additions
Eventually, this model will be expanded to add acidic calcium storages. To access the base model, make sure to select the correct branch. 

## Source
Bondarenko VE (2014) A Compartmentalized Mathematical Model of the β1-Adrenergic Signaling System in Mouse Ventricular Myocytes. PLoS ONE 9(2): e89113. https://doi.org/10.1371/journal.pone.0089113 
