==========================
Installation
==========================

EarthSHAB is a Python-based simulation framework for high-altitude solar balloon trajectory modeling.

**1. Clone the repository**

.. code-block:: bash

   git clone https://github.com/tkschuler/EarthSHAB.git
   cd EarthSHAB


**2. Create a virtual environment (optional but recommended)**

It's reccomended to install cfgrib with conda rather than pip

Using ``conda``:

.. code-block:: bash

   conda create -n earthshab python=3.11 pip
   conda activate earthshab
   conda install conda-forge::cfgrib

**3. Install dependencies**

.. code-block:: bash

   pip install -r requirements.txt
   pip install -e .


Platform Notes
--------------

- **Ubuntu / WSL (Recommended):**
  Fully supported and tested environment.

- **Windows:**
  Works best through WSL. Native Windows may require additional setup for some dependencies.

- **macOS:**
  Generally supported, but may require manual installation of system libraries for ``netCDF4``.

