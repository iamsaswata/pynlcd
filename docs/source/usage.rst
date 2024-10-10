Usage
=====

pynlcd can be used as a command-line tool or imported as a Python module.

Command-line Usage
------------------

To download NLCD data for a specific region:

.. code-block:: bash

   pynlcd --shapefile path/to/shapefile.shp --year 2021 --output path/to/output

For a specific point:

.. code-block:: bash

   pynlcd --point 37.7749 -122.4194 --year 2021 --output path/to/output

For a specific extent:

.. code-block:: bash

   pynlcd --extent -122.5 37.7 -122.3 37.8 --year 2021 --output path/to/output

Python Module Usage
-------------------

You can also use pynlcd in your Python scripts:

.. code-block:: python

   from pynlcd import get_land_cover
   from osgeo import ogr

   # Create your ROI geometry using GDAL/OGR
   shapefile_path = 'path/to/shapefile.shp'
   driver = ogr.GetDriverByName('ESRI Shapefile')
   dataSource = driver.Open(shapefile_path, 0)  # 0 means read-only. 1 means writeable.
   layer = dataSource.GetLayer()
   roi_geom = layer.GetNextFeature().GetGeometryRef()

   # Define the extent
   extent = (min_x, max_x, min_y, max_y)

   # Download the land cover data
   get_land_cover(roi_geom, extent, year=2021, spatial_resolution=0.0003, output_path='path/to/output')