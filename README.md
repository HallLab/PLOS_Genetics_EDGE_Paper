# PLOS_Genetics_EDGE_Paper
Scripts and results for simulated data

## Simulation Scripts

  * casecon-noise.pl for Case/Control simulation
  * quant-noise.pl for Quantitative outcome simulation
  * sim_updated.R - R script for simulation, used for the results found here.

## PBS Scripts

These are batch file scripts used to run on the compute cluster.

### runSimAll.pbs

Simulate case-control data using sim_updated.R data and run it in PLATO (using each of the 5 different encodings) on the fly, resulting in `data-all-1000`.

The nested structure of this folder is:

* MAF of 0.05, 0.1, 0.3, and 0.5
* Case-control ratio of 1, 2, and 3 combined with sample size of 2k, 10k, and 50k to get Case:Control simulations of:
  * 500:1500
  * 1000:1000
  * 1500:500
  * 2500:7500
  * 5000:5000
  * 5000:15000
  * 7500:2500
  * 10000:10000
  * 15000:5000
* PENDIFF values of 0.25, 0.33, and 0.45
* 57 Simulation Models
  * 6 main-effect models 6 encodings (REC, SUB, ADD, SUP, DOM, HET)
  * 21 interaction-only models based on all combinations of those 6 models (named, for example like `ADD_x_DOM`)
  * 21 main effects models based on all combinations of those 6 models (named, for example, like `ADD_+_DOM`)
  * 8 models specified by penetrance tables (HR-HR, HR-HET, HR-HA, HET-HET, HET-HA, XOR, Hyp, RHyp)
  * 1 NULL model
* PLATO results using 5 different encodings (add, dom, rec, codom, and weight (aka EDGE))
  * Each result compressed in its own folder
  * A text file reporting average pvalues across all simulations for each encoding

### RunSimSNR.pbs

Simulate case-control data using sim_updated.R and run it in PLATO (using each of 5 different encodings) on the fly, resulting in `data-SNR-1000`
Results in `data-SNR-1k-noMAF` are similar, but only at a MAF of 0.05.

The nested structure of this folder is:

* 5 different SNR values (0.01, 0.025, 0.05, 0.075, 0.10)
* 4 different MAF values(0.05, 0.10, 0.30, 0.50)
* 52 Simulation Models
  * Same as above, except 18 main effects models and 19 interaction models, skipping:
    * ADD_+_HET
    * ADD_+_SUB
    * DOM_+_DOM
    * SUB_x_SUB
    * SUP_x_HET
* PLATO results using 5 different encodings

### RunSim.pbs

Similar to RunSimAll.pbs, but using a single SNR of 0.05