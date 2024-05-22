<p align="center">
  <img src="docs/Knoten_Logo.svg" alt="Knoten" width=200> 
</p>

# Knoten

A library to leverage python wrapped Community Sensor Models (CSMs) for common spatial/sensor operations and testing.

## References:

- CSM (usgscsm): https://github.com/DOI-USGS/usgscsm
- Abstraction Layer for Ephemerides (ALE): https://github.com/DOI-USGS/ale
<hr>

## Overview

We currently use Knoten to help test our supported CSM implementations against well established [ISIS3 camera models](https://github.com/DOI-USGS/ISIS3). In short, The CSM standard, now at version 3.0.3, is a framework that provides a well-defined application program interface (API) for multiple types of sensors and has been widely adopted by remote sensing software systems (e.g. BAE's Socet GXP, Harris Corp.'s ENVI, Hexagon's ERDAS Imagine, and recently added to the NASA AMES Stereo Pipeline [ASP]). Our support for CSM is explained in this [abstract](https://www.hou.usra.edu/meetings/informatics2018/pdf/6040.pdf) and a recently submitted paper (not yet available). Currently, we support **Framing** and **Pushbroom** (line scanner) types of sensor models in the [usgscsm](https://github.com/DOI-USGS/usgscsm) library. 

A secondary requirement for our CSM implementation requires an ALE-generated Image Support Data (ISD). ISDs contain the [SPICE-derived](https://naif.jpl.nasa.gov/naif/toolkit.html) positional (and when needed velocity) description for each image. You can find several generated JSON-formatted examples [here](examples/data/)

Please see the **status report** below for the current instruments we have implemented and how well they match our ISIS3 camera models. In the near future, we will continue to address the pixel offsets we currently see. Both the CSM implementations (usgscsm) and ALE are currently in active development and both will be updated as needed to decrease these errors. Thus, none of the instruments have been tested enough for full production use.

<hr>

## Installing

You can install the latest build via conda. 

```
conda install -c usgs-astrogeology -c conda-forge knoten
```

You can also do a local install using the following steps within a clone of the repository.

1. Install the dependencies.  
   Note: creating the environment may take around and hour.  If conda is too slow, try using mamba instead.
    ```
    conda env create -f environment.yml
    ```

1. Install the package
    ```
    python setup.py install
    ```

<hr>

## Building Docs

To build the docs:
1. Open the docs directory: `cd docs`
1. Create (if not existing) and activate the environment (or use the main knoten environment if it is already created).
   
   ```sh
   conda env create -f environment.yml #If it does not exist already
   conda activate knoten-docs
   ```
1. Run `sphinx-build -b html . public` to build the docs in the `docs/public` directory (or change "public" to a subdirectory name of your choosing).
1. Browse to the `index.html` file in `docs/public` to open the docs.  If needed, run a utility like `http-server` in the `docs/public` directory to serve the docs on localhost.

<hr>
