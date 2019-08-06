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
|        MRO CTX        |      [link](examples/mro_ctx_isis_cmp.ipynb)      | No (gross errors) |     sample mean=0.8; line mean=-896     |     sample mean=0.8; line mean=-0.3     |
|        MEX HRSC       |      [link](examples/mex_hrsc_isis_cmp.ipynb)     |         No        |      sample mean=-80; line mean=-3      |       sample mean=80; line mean=3       |
|        LROC NAC       |      [link](examples/lrocnac_isis_cmp.ipynb)      | No (gross errors) |      sample mean=-40; line mean=-5      |                   n/a                   |
| Kaguya Terrain Camera |     [link](examples/kaguya_tc_isis_cmp.ipynb)     | No (gross errors) |     sample mean=-5; line mean=-1516     |      sample mean=5; line mean=1906      |
|   Messenger MDIS NAC  | [link](examples/messenger_mdisnac_isis_cmp.ipynb) |         No        |      sample mean=-0.2; line mean=-5     |       sample mean=0.3; line mean=5      |
|      Cassini ISS      |      [link](examples/cassini_isis_cmp.ipynb)      |         No        |    sample mean=-12.3; line mean=-0.1    |     sample mean=-0.2; line mean=0.9     |
|  Dawn Framing Camera  |      [link](examples/dawn_fc_isis_cmp.ipynb)      |     sub-pixel; in testing for production    | sample mean=-1.9e-03; line mean=1.8e-03 | sample mean=1.9e-03; line mean=-1.8e-03 |

The Difference column (CSM -> ISIS) represents the mean difference in pixels from running usgscsm's *image2ground* and then back to the camera using ISIS3's *campt* (ground2image). The Difference column (ISIS -> CSM) is simply the reverse starting with ISIS3 first with *campt* (image2ground) and then usgscsm's *ground2image*.