---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.11.1
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

<!-- #region -->
# Generating DRMs for GBM detectors

The main purpose for ```gbm_drm_gen``` is a fast and easy way to generate detector response matrices for the Fermi-GBM detectors in a pythonic way.




When we want to generate DRMS for GBM detectors, we simply need to obtain the proper data 
<!-- #endregion -->

```python
import matplotlib.pyplot as plt
import astropy.io.fits as fits

import numpy as np

%matplotlib inline

from jupyterthemes import jtplot

jtplot.style(context="talk", grid=False)

from gbm_drm_gen import DRMGen
from gbm_drm_gen.utils.package_data import get_path_of_data_file




```

## quick start


To create a DRM generator for TTE data, we need the TTE, CSPEC, and the TRIGDAT data files.
* The CSPEC data contains the output side of the DRM's energy bounds. 
* The TRIGDAT data contains the spacecraft orientation data.

```python
trigdat_file = get_path_of_data_file('example_data/glg_trigdat_all_bn110721200_v01.fit')
cspec_file = get_path_of_data_file('example_data/glg_cspec_n6_bn110721200_v00.pha')


# create the generator
gbm_n6_generator = DRMGen.from_128_bin_data(det_name= "n6",
                             time=0, # time relative to T0 or trigger time.
                             trigdat = trigdat_file,
                             mat_type = 2, # direct response + atmospheric scattering
                             cspecfile = cspec_file)
```

We can set the location of the source directly. The first run can be a bit slow as ```numba``` is used in the background to be very fast. 

```python
gbm_n6_generator.set_location(ra=329,dec = -38.2)
```

We can now checkout the matrix object created in the background:

```python
gbm_n6_generator.matrix
```

Or we can input RA and DEC to create a 3ML style OGIP response directly:

```python
response = gbm_n6_generator.to_3ML_response(ra=329,dec = -38.2)
```

```python
fig = response.plot_matrix()
```

To see how the effective area varies with location, we can loop through various angles.

```python
fig, ax = plt.subplots()

bounds = np.vstack((gbm_n6_generator.monte_carlo_energies[:-1],gbm_n6_generator.monte_carlo_energies[1:])).T
de = np.diff(bounds)
ene = np.mean(bounds,axis=1)

for ra in np.linspace(260, 350, 10):
    
    gbm_n6_generator.set_location(ra=ra,dec = -38.2)

    ax.loglog(ene,gbm_n6_generator.matrix.sum(axis=0),label=r'%d$^{\circ}$'%ra)

ax.set_ylim(1)
ax.legend()
ax.set_xlabel(r'Effective Area (cm$^2$)');
ax.set_xlabel('MC Energies');
```

<!-- #region -->
## Into the details

Ok, now let's go through the various specifics of the DRM generator constructor.


First, the generator needs to know:

* The current location and oreintation of GBM for the data of interest.

* The channel to PHA reconstructed energy from the CSPEC files. 


For the first, one can get the current spacecraft from the either the triggers trigdat file, or for non-trggered data, one can obtain the position history file. These are available at the [NASA database](https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/). If using a position history file. You need to specify ```T0=<Fermi MET>``` in the constructor so that the time coordinate will be relative to this MET.

Internally, the class uses [gbmgeometry](https://gbmgeometry.readthedocs.io/en/latest/) to convert RA,Dec to the approriate spacecraft coordinates. However, one can also create DRMs in spacecraft coordinates directly. 

<!-- #endregion -->

## custom energy binning

Maybe you are a curious person and want to investigate a response with finer input energies to model line features in solar flares? 

It is possible to add a custom array of input energies. To do this we need to import the ```NaITTEEdges``` and ```BgoTTEEdges``` classes.


```python
from gbm_drm_gen import  NaiTTEEdges, BgoTTEEdges
```

These objects allow you to specify input (monte carlo) energies either from an array or in log spaced binning.

<div class="alert alert-info">

**Note:** The number of energies must be off and include the low and high end points which are specific to the NaI [5.0, 50000.0] and BGO [100., 200000.0] detectors.

</div>

### Log-spaced energies

First we will try with **VERY** fine log spaced binning above the typical 140 input energies


```python
custom_edges = NaiTTEEdges.from_log_bins(n_bins=531)


gbm_n6_generator = DRMGen.from_128_bin_data(det_name= "n6",
                             trigdat = trigdat_file,
                             mat_type = 2, # direct response + atmospheric scattering
                             cspecfile = cspec_file,
                             # pass the custom edges
                             custom_input_edges=custom_edges
                            
                            )


response = gbm_n6_generator.to_3ML_response(ra=329,dec = -38.2)

fig = response.plot_matrix()
```

And now with a much coarser input binning:

```python
custom_edges = NaiTTEEdges.from_log_bins(n_bins=73)


gbm_n6_generator =  DRMGen.from_128_bin_data(det_name= "n6",
                             trigdat = trigdat_file,
                             mat_type = 2, # direct response + atmospheric scattering
                             cspecfile = cspec_file,
                             # pass the custom edges
                             custom_input_edges=custom_edges
                            
                            )


response = gbm_n6_generator.to_3ML_response(ra=329,dec = -38.2)

fig = response.plot_matrix()
```

It is easier to see the difference with simple matrix plotting:

```python
custom_edges = NaiTTEEdges.from_log_bins(n_bins=91)


gbm_n6_generator =  DRMGen.from_128_bin_data(det_name= "n6",
                             trigdat = trigdat_file,
                             mat_type = 2, # direct response + atmospheric scattering
                             cspecfile = cspec_file,
                             # pass the custom edges
                             custom_input_edges=custom_edges
                            
                            )

gbm_n6_generator.set_location(ra=329,dec = -38.2)


plt.matshow(gbm_n6_generator.matrix.T)

custom_edges = NaiTTEEdges.from_log_bins(n_bins=141)


gbm_n6_generator =  DRMGen.from_128_bin_data(det_name= "n6",
                             trigdat = trigdat_file,
                             mat_type = 2, # direct response + atmospheric scattering
                             cspecfile = cspec_file,
                             # pass the custom edges
                             custom_input_edges=custom_edges
                            
                            )

gbm_n6_generator.set_location(ra=329,dec = -38.2)


plt.matshow(gbm_n6_generator.matrix.T)



custom_edges = NaiTTEEdges.from_log_bins(n_bins=1541)


gbm_n6_generator =  DRMGen.from_128_bin_data(det_name= "n6",
                             trigdat = trigdat_file,
                             mat_type = 2, # direct response + atmospheric scattering
                             cspecfile = cspec_file,
                             # pass the custom edges
                             custom_input_edges=custom_edges
                            
                            )

gbm_n6_generator.set_location(ra=329,dec = -38.2)


plt.matshow(gbm_n6_generator.matrix.T)
```

### a custom array

And we can even supply and entirely custom array of energies:

```python
edges = np.linspace(5., 50000., 1001)


custom_edges = NaiTTEEdges.from_custom_array(edges)


gbm_n6_generator =  DRMGen.from_128_bin_data(det_name= "n6",
                             trigdat = trigdat_file,
                             mat_type = 2, # direct response + atmospheric scattering
                             cspecfile = cspec_file,
                             # pass the custom edges
                             custom_input_edges=custom_edges
                            
                            )


response = gbm_n6_generator.to_3ML_response(ra=329,dec = -38.2)

fig = response.plot_matrix()
```

## Creating RSP2 files

```python
from gbm_drm_gen import create_rsp2
```

```python
# create the generator
gbm_n6_generator = DRMGen.from_128_bin_data(det_name= "n6",
                             time=0, # time relative to T0 or trigger time.
                             trigdat = trigdat_file,
                             mat_type = 2, # direct response + atmospheric scattering
                             cspecfile = cspec_file)
```

```python
output_file_name = "my_new_rsp.rsp2" # you must call it an RSP2 file!
```

```python
create_rsp2(output_file_name,
            response_generator=gbm_n6_generator,
            ra=0,
            dec=0,
            tstart=0,
            tstop=10,
            delta_time=2,
            overwrite=True
           )
```

```python
with fits.open(output_file_name) as f:
    
    f.info()
    
    print(f[0].header['DRM_NUM'])
    print(f[2].header['TSTART'])
    print(f[2].header['TSTOP'])
    
   
```

```python
from threeML import TimeSeriesBuilder
```

```python
tsb = TimeSeriesBuilder.from_gbm_cspec_or_ctime("cspec", cspec_or_ctime_file=cspec_file, rsp_file="my_new_rsp.rsp2")

tsb.view_lightcurve();

```

```python
tsb.set_active_time_interval('1-5')
```

```python

```
