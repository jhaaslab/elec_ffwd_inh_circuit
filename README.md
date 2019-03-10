# Electrical synapses and transient signals in feedforward canonical circuits  
This is the code to run simulations in the paper:  
*Full title*: **Electrical synapses regulate both subthreshold and population activity of principal cells in response to transient inputs within canonical feedforward circuits**  
Link to paper on *PLoS Comput Biol*: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006440   
Prepint of submission on *bioRxiv*: https://www.biorxiv.org/content/early/2018/08/17/394593   
*Authors*: Tuan Pham (1) and Julie S. Haas (2)   
*Affiliation*: Dept. of Biol. Sci., Lehigh University, Bethlehem, PA, USA  
*Contact*: (1) phhanhtuan*AT*gmail*DOT*com, (2) julie*DOT*haas*AT*lehigh*DOT*edu, jhaaslab*AT*gmail*DOT*com  

### Requirements
+ `MATLAB R2018a` 
+ `Python 2.7` 
+ [`Brian 2`](https://brian2.readthedocs.io/en/stable/)

### File organization
+ `simulations` contains simulation scripts in the following `jupyter` notebooks:  
	+ Subthreshold integration in canonical circuits: `subthreshold_simumations.ipynb` 
	+ Spiking responses in networks of canonical circuits: `network_activity_simulations.ipynb`
+ `prm_files` contains neuron model parameters and network configuration parameters
+ `plotting_figures` contains the `MATLAB` scripts for analysis and figure generation, using the functions in `functions` folder
+ `figures` contains generated figures in the paper 
