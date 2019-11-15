# Knoten

A library to leverage python wrapped Community Sensor Models (CSMs) for common spatial/sensor operations and testing.

**References:**

- CSM (usgscsm): https://github.com/USGS-Astrogeology/usgscsm
- Abstraction Layer for Ephemerides (ALE): https://github.com/USGS-Astrogeology/ale
<hr>

**Overview:**

We currently use Knoten to help test our supported CSM implementations against well established [ISIS3 camera models](https://github.com/USGS-Astrogeology/ISIS3). In short, The CSM standard, now at version 3.0.3, is a framework that provides a well-defined application program interface (API) for multiple types of sensors and has been widely adopted by remote sensing software systems (e.g. BAE's Socet GXP, Harris Corp.'s ENVI, Hexagon's ERDAS Imagine, and recently added to the NASA AMES Stereo Pipeline [ASP]). Our support for CSM is explained in this [abstract](https://www.hou.usra.edu/meetings/informatics2018/pdf/6040.pdf) and a recently submitted paper (not yet available). Currently, we support **Framing** and **Pushbroom** (line scanner) types of sensor models in the [usgscsm](https://github.com/USGS-Astrogeology/usgscsm) library. 

A secondary requirement for our CSM implementation requires an ALE-generated Image Support Data (ISD). ISDs contain the [SPICE-derived](https://naif.jpl.nasa.gov/naif/toolkit.html) positional (and when needed velocity) description for each image. You can find several generated JSON-formatted examples [here](examples/data/)

Please see the **status report** below for the current instruments we have implemented and how well they match our ISIS3 camera models. In the near future, we will continue to address the pixel offsets we currently see. Both the CSM implementations (usgscsm) and ALE are currently in active development and both will be updated as needed to decrease these errors. Thus, none of the instruments have been tested enough for full production use.

<hr>

Installing:

1. download knoten's environment.yml
2. Install Anaconda or Miniconda
3. Within an conda terminal, type:
```
conda env create -f environment.yml
```

<hr>

**Status Report - July 2019:**

For full testing reports and example usage, please see the linked example Jupyter notebooks in the table below. 

|       Instrument      |                      Jupyter Notebooks                     |  Production Ready |    Difference CSM -> ISIS (in pixels)   |    Difference ISIS -> CSM (in pixels)   |
|:---------------------:|:-------------------------------------------------:|:-----------------:|:---------------------------------------:|:---------------------------------------:|
|       MRO HiRISE      |     [link](examples/mro_hirise_isis_cmp.ipynb)    |     sub-pixel; in testing for production    | sample mean=-2.0e-05; line mean=2.5e-08 | sample mean=-3.0e-08; line mean=1.2e-04 |
|        MRO CTX        |      [link](examples/mro_ctx_isis_cmp.ipynb)      | nearly sub-pixel; in testing |          |     sample mean=1.4; line mean=0.01     |
|        MEX HRSC       |      [link](examples/mex_hrsc_isis_cmp.ipynb)     |         sub-pixel; in testing for production        |      sample mean=-80; line mean=-3      |       sample mean=0.0004 ; line mean=0.05       |
|        LROC NAC       |      [link](examples/lrocnac_isis_cmp.ipynb)      | sub-pixel; in testing for production |      sample mean=-0.003; line mean=-0.0006      |                   sample mean=0.003	line mean=0.0005                   |
| Kaguya Terrain Camera |     [link](examples/kaguya_tc_isis_cmp.ipynb)     | barely sub-pixel; in testing |     sample mean=0.5; line mean=-0.9     |      sample mean=-0.5; line mean=0.94      |
|   Messenger MDIS NAC  | [link](examples/messenger_mdisnac_isis_cmp.ipynb) |  sub-pixel; in testing for production |      sample mean=-0.01; line mean=-0.003     |       sample mean=0.01; line mean=0.003      |
|  Cassini ISS NAC      |      [link](examples/cassini_isis_nac_cmp.ipynb)      |  sub-pixel; in testing for production  |    sample mean=-0.001; line mean=0.01    |     sample mean=0.001; line mean=-0.01     |
|  Cassini ISS WAC      |      [link](examples/cassini_isis_wac_cmp.ipynb)      |  sub-pixel; in testing for production  |    sample mean=0.001; line mean=0.004    |     sample mean=-0.001; line mean=-0.004     |
|  Dawn Framing Camera  |      [link](examples/dawn_fc_isis_cmp.ipynb)      |     sub-pixel; in testing for production    | sample mean=-0.02; line mean=0.003 | sample mean=0.02; line mean=-0.003 |

The Difference column (CSM -> ISIS) represents the mean difference in pixels from running usgscsm's *image2ground* and then back to the camera using ISIS3's *campt* (ground2image). The Difference column (ISIS -> CSM) is simply the reverse starting with ISIS3 first with *campt* (image2ground) and then usgscsm's *ground2image*.
