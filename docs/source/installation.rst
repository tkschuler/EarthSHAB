==========================
Installation
==========================

EarthSHAB is a Python-based simulation framework for high-altitude solar balloon trajectory modeling.

**1. Clone the repository**

.. code-block:: bash

   git clone https://github.com/tkschuler/EarthSHAB.git
   cd EarthSHAB


**2. Create a virtual environment (optional but recommended)**

Using ``conda``:

.. code-block:: bash

   conda create -n earthshab python=3.11
   conda activate earthshab
   pip install -r requirements.txt
   pip install -e .

**3. Install dependencies**

.. code-block:: bash

   pip install -r requirements.txt


**4. Install EarthSHAB (editable mode)**

.. code-block:: bash

   pip install -e .


Platform Notes
--------------

- **Ubuntu / WSL (Recommended):**
  Fully supported and tested environment.

- **Windows:**
  Works best through WSL. Native Windows may require additional setup for some dependencies.

- **macOS:**
  Generally supported, but may require manual installation of system libraries for ``netCDF4``.

