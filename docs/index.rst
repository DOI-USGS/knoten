Knoten
======

A library to leverage python wrapped Community Sensor Models (CSMs) for common spatial/sensor operations and testing.


Table of Contents
-----------------

:Date: |today|

.. toctree::
   :maxdepth: 2

   library/index
   

References
----------

- CSM (usgscsm): https://github.com/DOI-USGS/usgscsm
- Abstraction Layer for Ephemerides (ALE): https://github.com/DOI-USGS/ale

Overview
--------

We currently use Knoten to help test our supported CSM implementations against well established `ISIS3 camera models <https://github.com/DOI-USGS/ISIS3>`_. In short, The CSM standard, now at version 3.0.3, is a framework that provides a well-defined application program interface (API) for multiple types of sensors and has been widely adopted by remote sensing software systems (e.g. BAE's Socet GXP, Harris Corp.'s ENVI, Hexagon's ERDAS Imagine, and recently added to the NASA AMES Stereo Pipeline [ASP]). Our support for CSM is explained in this `abstract <https://www.hou.usra.edu/meetings/informatics2018/pdf/6040.pdf>`_ and a recently submitted paper (not yet available). Currently, we support **Framing** and **Pushbroom** (line scanner) types of sensor models in the `usgscsm <https://github.com/DOI-USGS/usgscsm>`_ library. 

A secondary requirement for our CSM implementation requires an ALE-generated Image Support Data (ISD). ISDs contain the `SPICE-derived <https://naif.jpl.nasa.gov/naif/toolkit.html>`_ positional (and when needed velocity) description for each image. You can find several generated JSON-formatted examples `here <https://github.com/DOI-USGS/knoten/blob/main/examples/data>`_.

Please see the **status report** below for the current instruments we have implemented and how well they match our ISIS3 camera models. In the near future, we will continue to address the pixel offsets we currently see. Both the CSM implementations (usgscsm) and ALE are currently in active development and both will be updated as needed to decrease these errors. Thus, none of the instruments have been tested enough for full production use.


Installing
----------

You can install the latest build via conda::

   conda install -c usgs-astrogeology -c conda-forge knoten

You can also do a local install using the following steps within a clone of the repository.

#. Install the dependencies.  
   Note: creating the environment may take around and hour.  If conda is too slow, try using mamba instead::
   
      conda env create -f environment.yml

#. Install the package::

      python setup.py install

Building Docs
-------------

To build the docs:

#. Open the docs directory: ``cd docs``
#. Create (if not existing) and activate the environment (or use the main knoten environment if it is already created)::

      conda env create -f environment.yml #If it does not exist already
      conda activate knoten-docs

#. Run ``sphinx-build -b html . public`` to build the docs in the ``docs/public`` directory (or change "public" to a subdirectory name of your choosing).
#. Browse to the ``index.html`` file in ``docs/public`` to open the docs.  If needed, run a utility like ``http-server`` in the ``docs/public`` directory to serve the docs on localhost.


Status Report - November 2019
-----------------------------

For full testing reports and example usage, please see the linked example Jupyter notebooks in the table below. 

 ======================= ============================================================= ====================================== =========================================== ============================================== 
  Instrument              Jupyter Notebooks                                             Production Ready                       Difference CSM -> ISIS (in pixels)          Difference ISIS -> CSM (in pixels)            
 ======================= ============================================================= ====================================== =========================================== ============================================== 
  MRO HiRISE              `MRO HiRISE <examples/mro_hirise_isis_cmp.ipynb>`_            sub-pixel; in testing for production   sample mean=-2.0e-05; line mean=2.5e-08     sample mean=-3.0e-08; line mean=1.2e-04       
  MRO CTX                 `MRO CTX <examples/mro_ctx_isis_cmp.ipynb>`_                  nearly sub-pixel; still in research    gross error in line                         sample mean=0.0002; line mean=-0.07           
  MEX HRSC                `MEX HRSC <examples/mex_hrsc_isis_cmp.ipynb>`_                sub-pixel; in testing for production   sample mean=0.000038; line mean=-0.000072   sample mean=-0.000038 ; line mean=-7.512e-05  
  LROC NAC                `LROC NAC <examples/lrocnac_isis_cmp.ipynb>`_                 sub-pixel; in testing for production   sample mean=-0.003; line mean=-0.0006       sample mean=0.0005	line mean=0.003          
  Kaguya Terrain Camera   `Kaguya TC <examples/kaguya_tc_isis_cmp.ipynb>`_              barely sub-pixel; in testing           sample mean=0.0001; line mean=0.00003       sample mean=0.009; line mean=-1.242           
  Messenger MDIS NAC      `MGR MDIS NAC <examples/messenger_mdisnac_isis_cmp.ipynb>`_   sub-pixel; in testing for production   sample mean=-0.01; line mean=-0.003         sample mean=0.01; line mean=0.003             
  Cassini ISS NAC         `CAS ISS NAC <examples/cassini_isis_nac_cmp.ipynb>`_          sub-pixel; in testing for production   sample mean=-0.001; line mean=0.01          sample mean=0.001; line mean=-0.01            
  Cassini ISS WAC         `CAS ISS WAC <examples/cassini_isis_wac_cmp.ipynb>`_          sub-pixel; in testing for production   sample mean=0.001; line mean=0.004          sample mean=-0.001; line mean=-0.004          
  Dawn Framing Camera     `Dawn FC <examples/dawn_fc_isis_cmp.ipynb>`_                  sub-pixel; in testing for production   sample mean=-0.02; line mean=0.003          sample mean=0.02; line mean=-0.003            
 ======================= ============================================================= ====================================== =========================================== ==============================================  

The Difference column (CSM -> ISIS) represents the mean difference in pixels from running usgscsm's *image2ground* and then back to the camera using ISIS3's *campt* (ground2image). The Difference column (ISIS -> CSM) is simply the reverse starting with ISIS3 first with *campt* (image2ground) and then usgscsm's *ground2image*.
