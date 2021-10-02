<p align="center">
    <img alt="Illustration of SPOD workflow" src="docs/figs/logo.svg" width="400" />
</p>

---

# Turbulence Analyzer in Python

This repository contains a Python toolkit that calculates and visualizes turbulence anisotropy and turbulent viscosity from Reynolds stress components. A more detailed explanation can be found in the [theory guide](./docs/theory_guide.pdf). A graphical abstract of this package is illustrated in the following figure.

<p align="center">
    <img alt="Illustration of TurbAna" src="docs/figs/TurbAna_schematic.png" width="500" />
</p>

If this script appears useful for your research, an explicit mention of the work [[5](#ddes-comp)] (for turbulence anisotropy) and [[6](#ddes-comp)] (for turbulent viscosity) would be highly appreciated.

## Quick Start


### Step 1: Download package
#### Download from Git clone in the terminal

git clone [https://github.com/HexFluid/TurbAna.git](https://github.com/HexFluid/TurbAna.git)

#### Download from browser

Download from this [link](https://github.com/HexFluid/TurbAna/archive/master.zip) and then unzip it.

### Step 2: Install prerequisites
Launch a terminal (UNIX) or an Anaconda Prompt (Windows) window and change directory to *TurbAna*. Run the following command line to install/upgrade the prerequisite Python packages.

```
pip install -r requirements.txt
```

### Step 3: Load example data
Run the following script with Python 3 to load the data:
```python
import h5py
import os

current_path = os.getcwd() # assuming Python launched in the 'TurbAna' dir
data_path    = os.path.join(current_path,'tutorials','bstep_data','bstep_DDES.h5')

h5f  = h5py.File(data_path,'r')
data = h5f['data'][:]        # flow field data
h5f.close()
```

### Step 4: Calculate and visualize turbulence anisotropy
Run the following script to obtain SPOD results:
```python
import spod

spod.spod(data,dt,current_path,weight='default',nOvlp='default',window='default',method='fast')

SPOD_LPf  = h5py.File(os.path.join(current_path,'SPOD_LPf.h5'),'r') # load data from h5 format
L = SPOD_LPf['L'][:,:]    # modal energy E(f, M)
P = SPOD_LPf['P'][:,:,:]  # mode shape
f = SPOD_LPf['f'][:]      # frequency
SPOD_LPf.close()
```

### Step 5: Calculate and visualize turbulence anisotropy
Finally, run the following script to visualize the SPOD spectrum:
```python
fig = spod.plot_spectrum(f,L,hl_idx=5)
```

Expected results:
<p align="left">
    <img alt="Illustration of SPOD workflow" src="docs/figs/SPOD_quickstart_result.png" width="400" />
</p>

For more postprocess tutorials including plotting mode shapes and reconstructed flow fields, please refer to the scripts with detailed comments in [tutorials](./tutorials/README.md).



## List of Files

<pre>
.
|-- docs
|   |-- figs
|   |-- theory_guide.pdf
|-- tutorials
|   |-- bstep_data
|   |   |-- results
|   |   |-- bstep_DDES.h5
|   |   |-- bstep_DNS.h5
|   |   |-- bstep_SA_frozen.h5
|   |-- bump_data
|   |   |-- results
|   |   |-- bump_LES.h5
|   |   |-- bump_SST.h5
|   |-- cooling_data
|   |   |-- results
|   |   |-- cooling_DDES.h5
|   |-- SBLI_data
|   |   |-- results
|   |   |-- SBLI_DNS.h5
|   |   |-- SBLI_SST.h5
|   |-- 01_bstep.py
|   |-- 02_bump.py
|   |-- 03_cooling.py
|   |-- 04_SBLI.py
|-- LICENSE
|-- requirements.txt
|-- TurbAna.py
</pre>

- **TurbAna.py**: main script of turbulence analyzer
- **requirements.txt**: a list of prerequisite Python libraries
- **LICENSE**: license file
- **tutorials**
  - **bstep_data**
    - **results**: postprocess results of the single variable case
    - **bstep_DDES.h5**: flow field and grid data (HDF5 format)
    - **bstep_DNS.h5**: flow field and grid data (HDF5 format)
    - **bstep_SA_frozen.h5**: flow field and grid data (HDF5 format)
    - **bstep_SAQCR_frozen.h5**: flow field and grid data (HDF5 format)
    - **bump_data**
      - **results**: postprocess results of the single variable case
      - **bump_LES.h5**: flow field and grid data (HDF5 format)
      - **bump_SST.h5**: flow field and grid data (HDF5 format)
  - **cooling_data**
    - **results**: postprocess results of the multiple variable case
    - **cooling_DDES.h5**: flow field and grid data (HDF5 format)
  - **SBLI_data**
    - **results**: postprocess results of the single variable case
    - **SBLI_DNS.h5**: flow field and grid data (HDF5 format)
    - **SBLI_SST.h5**: flow field and grid data (HDF5 format)
  - **01_bstep.py**: tutorial script for the backward-facing step case
  - **02_bump.py**: tutorial script for the transonic bump case
  - **03_cooling.py**: tutorial script for the film cooling case
  - **04_SBLI.py**: tutorial script for the shock boundary layer interaction (SBLI) case
- **docs**
  - **figs**: figures appeared in the markdown files
  - **theory_guide.pdf**: theory of TurbAna

## References
[<a id="anisotropy">1</a>] Emory, M., & Iaccarino, G. (2014). Visualizing turbulence anisotropy in the spatial domain with componentality contours. Center for Turbulence Research Annual Research Briefs, 123-138. [[link](https://web.stanford.edu/group/ctr/ResBriefs/2014/14_emory.pdf)]

[<a id="ddes-comp">2</a>] He, X., Zhao, F., & Vahdati, M. (2022). Detached eddy simulation: recent development and application to compressor tip leakage flow. ASME Journal of Turbomachinery, 144(1), 011009. [[DOI](https://doi.org/10.1115/1.4052019)][[preprint](https://www.researchgate.net/publication/347355348_Detached_Eddy_Simulation_Recent_Development_and_Application_to_Compressor_Tip_Leakage_Flow)]

[<a id="SA-ML">3</a>] He, X., Tan, J., & Vahdati, M. (2021). Towards Explainable Machine Learning Assisted Turbulence Modelling for Transonic Flows. [[preprint](https://www.researchgate.net/publication/344903748_Towards_Explainable_Machine_Learning_Assisted_Turbulence_Modelling_for_Transonic_Flows)]
