from osgeo import gdal
from osgeo import ogr
from osgeo import osr
import os
import math
import sys
import argparse
import multiprocessing
from functools import partial

def main():
  parser = argparse.ArgumentParser(description='Download NLCD data for specified year and region.')
  group = parser.add_mutually_exclusive_group(required=True)
  group.add_argument('--shapefile', help='Path to the input shapefile defining the region of interest')
  group.add_argument('--point', nargs=2, metavar=('LAT', 'LON'), type=float, help='Point coordinates (latitude longitude)')
  group.add_argument('--extent', nargs=4, metavar=('MINX', 'MINY', 'MAXX', 'MAXY'), type=float, help='Extent defined by minx miny maxx maxy')

  parser.add_argument('--year', type=int, default=2021, help='Year of NLCD data to download (default: 2021)')
  parser.add_argument('--cellsize', type=float, default=None, help='Cell size in degrees (optional)')
  parser.add_argument('--output', required=True, help='Output path where tiles and combined tif will be saved')
  args = parser.parse_args()

  # Validate the year
  valid_years = [2001, 2004, 2006, 2008, 2011, 2013, 2016, 2019, 2021]
  if args.year not in valid_years:
      print(f"Error: Year {args.year} is not valid. Valid years are: {valid_years}")
      sys.exit(1)

  # Get ROI geometry
  roi_geom = get_roi_geometry(args)
  # Get the bounding box of the ROI
  min_x, max_x, min_y, max_y = roi_geom.GetEnvelope()
  extent = (min_x, max_x, min_y, max_y)

  # Create output directory if it doesn't exist
  if not os.path.exists(args.output):
      os.makedirs(args.output)

  get_land_cover(roi_geom, extent, args.year, args.cellsize, args.output)

  # Crop the TIFF using the shapefile if provided
  if args.shapefile:
    crop_tif_with_shapefile(args.output, args.shapefile, args.year)  


  # Handle point input
  if args.point:
    create_point_file(args.output, args.point, args.year)
    return    

def create_point_file(output_path, point, year):
    """
    Create a file containing the NLCD value for a specific point and year.
    This function retrieves the NLCD (National Land Cover Database) value for a given latitude and longitude
    at a specified year, writes this information to a text file, and then removes the corresponding TIFF file.

    :param output_path: str
        The directory where the output text file and TIFF file are located.

    :param point: tuple
        A tuple containing the latitude and longitude of the point (lat, lon).

    :param year: int
        The year for which the NLCD value is to be retrieved.

        
    :return: None

    """

    lat, lon = point
    nlcd_value = get_nlcd_value_at_point(output_path, lat, lon, year)
    point_file = os.path.join(output_path, f"NLCD_{year}_{lat}_{lon}.txt")
    with open(point_file, 'w') as f:
        f.write(f"Latitude: {lat}\nLongitude: {lon}\nNLCD Value: {nlcd_value}\n")
    print(f"Point coordinates and NLCD value saved to {point_file}")

    # remove the tif file
    tif_file = os.path.join(output_path, f"NLCD_{year}_Land_Cover.tif")
    os.remove(tif_file)

def get_nlcd_value_at_point(output_path, lat, lon, year):
    """
    Retrieve the NLCD (National Land Cover Database) value at a specific geographic point.

    Parameters
    ----------

    output_path : str
        The directory path where the NLCD TIFF files are stored.
    lat : float
        The latitude of the point of interest.
    lon : float
        The longitude of the point of interest.
    year : int
        The year of the NLCD data to be used.

    Returns
    -------
    nlcd : int
        The NLCD value at the specified geographic point.
    """
    
    # Path to the combined TIFF file
    combined_tif = os.path.join(output_path, f"NLCD_{year}_Land_Cover.tif")

    # Open the TIFF file
    dataset = gdal.Open(combined_tif)
    if not dataset:
        print(f"Error: Unable to open {combined_tif}")
        sys.exit(1)

    # Get the geotransform and inverse geotransform
    transform = dataset.GetGeoTransform()
    inv_transform = gdal.InvGeoTransform(transform)

    # Convert lat/lon to pixel coordinates
    px, py = gdal.ApplyGeoTransform(inv_transform, lon, lat)

    # Read the pixel value
    band = dataset.GetRasterBand(1)
    nlcd_value = band.ReadAsArray(int(px), int(py), 1, 1)[0][0]

    return nlcd_value



def crop_tif_with_shapefile(output_path, shapefile_path, year):
    """
    Crop the combined TIFF file using the provided shapefile.
    This function uses the provided shapefile to crop the combined TIFF file
    and saves the cropped TIFF file to the output directory.

    Parameters
    ----------
    output_path : str
        The directory path where the combined and cropped TIFF files are stored.
    shapefile_path : str
        The path to the shapefile used for cropping.
    year : int
        The year of the NLCD data to be cropped.

    Returns
    -------
    None
    """

    # Path to the combined TIFF file
    combined_tif = os.path.join(output_path, f"NLCD_{year}_Land_Cover.tif")
    # Path to the cropped TIFF file
    cropped_tif = os.path.join(output_path, f"NLCD_{year}_Land_Cover_final.tif")

    # Use gdal.Warp to crop the TIFF using the shapefile
    gdal.Warp(cropped_tif, combined_tif, cutlineDSName=shapefile_path, cropToCutline=True)

    # Remove the combined TIFF file
    os.remove(combined_tif)

    # Rename the cropped TIFF file to the final name
    os.rename(cropped_tif, combined_tif)

    print(f"Cropped TIFF saved to {cropped_tif}")    


def get_land_cover(roi_geom, extent, year, spatial_resolution, output_path):
  """
  Download NLCD data for the specified year and region.

  Parameters
  ----------
  roi_geom : ogr.Geometry
    An OGR geometry representing the region of interest.

  extent : tuple
    A tuple containing the minimum and maximum x and y coordinates of the extent.

  year : int
    The year of the NLCD data to be downloaded.

  spatial_resolution : float
    The spatial resolution of the data in degrees.

  output_path : str
    The directory where the downloaded data will be saved.

  Returns
  -------
  None
  """

  
  min_x, max_x, min_y, max_y = extent
  span_x = abs(max_x - min_x)
  span_y = abs(max_y - min_y)
  area = span_x * span_y

  # Decide whether to download in tiles or as a whole
  if area <= 9:  # 9 square degrees
      # Just get the whole image
      out_img = os.path.join(output_path, f"NLCD_{year}_Land_Cover.tif")
      wcs_url = get_url(extent, year)
      print("Getting the image...")
      result = get_IMG(wcs_url, out_img, spatial_resolution)
  else:
      # Split extent into tiles of ~1 sq deg
      tile_extents = split_extent(extent, spatial_resolution, roi_geom)
      if not tile_extents:
          print("No tiles intersect with the ROI.")
          return
      failed_tiles = []
      tile_paths = []

      # Use multiprocessing Pool for parallel downloads
      num_processes = min(len(tile_extents), multiprocessing.cpu_count())
      pool = multiprocessing.Pool(processes=num_processes)
      print(f"Downloading {len(tile_extents)} tiles using {num_processes} processes...")

      # Prepare arguments for parallel processing
      download_args = []
      for i, tile_ext in enumerate(tile_extents):
          out_tile = os.path.join(output_path, f"tile_{i}.tif")
          wcs_url = get_url(tile_ext, year)
          download_args.append((wcs_url, out_tile, spatial_resolution))

      # Download tiles in parallel
      results = pool.map(download_worker, download_args)
      pool.close()
      pool.join()

      # Collect successful tile paths
      for success, out_tile in results:
          if success:
              tile_paths.append(out_tile)
          else:
              failed_tiles.append(out_tile)

      if failed_tiles:
          print("Retrying failed tiles...")
          # Retry failed downloads
          pool = multiprocessing.Pool(processes=num_processes)
          retry_args = [(get_url(tile_extents[i], year), failed_tiles[i], spatial_resolution)
                        for i in range(len(failed_tiles))]
          retry_results = pool.map(download_worker, retry_args)
          pool.close()
          pool.join()

          for success, out_tile in retry_results:
              if success:
                  tile_paths.append(out_tile)
              else:
                  print(f"Failed to download tile: {out_tile}")

      # Mosaic tiles together
      if tile_paths:
          output_mosaic = os.path.join(output_path, f"NLCD_{year}_Land_Cover.tif")
          build_mosaic(tile_paths, output_mosaic)
      else:
          print("No tiles were successfully downloaded.")


def download_worker(args):
  """
  Download a single tile using the given WCS URL and save it to the output path.

  Parameters
  ----------
  args : tuple
    A tuple containing the WCS URL, output image path, and spatial resolution.

  Returns
  -------
  tuple
    A tuple containing a boolean indicating success or failure and the output image path.

 """
  
  wcs_url, out_img, spatial_resolution = args
  result = get_IMG(wcs_url, out_img, spatial_resolution)
  return (result == "success", out_img)


def split_extent(extent, spatial_resolution, roi_geom):
  """
  
  Split the extent into subextents of approximately 1 square degree each.
  This function divides the extent into subextents of approximately 1 square degree each,
  taking into account the spatial resolution and the overlap between tiles.

  Parameters
  ----------
  extent : tuple
    A tuple containing the minimum and maximum x and y coordinates of the extent.

  spatial_resolution : float
    The spatial resolution of the data in degrees.

  roi_geom : ogr.Geometry
    An OGR geometry representing the region of interest.

  Returns
  -------
  list
    A list of subextents that intersect with the region of interest

  """
  # Size of each subextent (approximately 1 square degree)
  min_x, max_x, min_y, max_y = extent
  subextent_size = 1.0  # in degrees
  # Add some tile overlap to prevent possible gaps
  if spatial_resolution:
      tile_margin = spatial_resolution
  else:
      tile_margin = 0.0006  # 0.0006 degrees = 60 meters

  subextents = []

  # Calculate the number of rows and columns
  num_rows = math.ceil((max_y - min_y) / subextent_size)
  num_cols = math.ceil((max_x - min_x) / subextent_size)

  step_x = (max_x - min_x) / num_cols
  step_y = (max_y - min_y) / num_rows

  # Create spatial reference for 4326
  srs = osr.SpatialReference()
  srs.ImportFromEPSG(4326)

  for row in range(num_rows):
      for col in range(num_cols):
          # Calculate the bounds of the subextent
          subextent_min_x = (min_x + col * step_x) - tile_margin
          subextent_max_x = (min_x + (col + 1) * step_x) + tile_margin
          subextent_min_y = (min_y + row * step_y) - tile_margin
          subextent_max_y = (min_y + (row + 1) * step_y) + tile_margin

          # Create tile geometry
          ring = ogr.Geometry(ogr.wkbLinearRing)
          ring.AddPoint(subextent_min_x, subextent_min_y)
          ring.AddPoint(subextent_min_x, subextent_max_y)
          ring.AddPoint(subextent_max_x, subextent_max_y)
          ring.AddPoint(subextent_max_x, subextent_min_y)
          ring.AddPoint(subextent_min_x, subextent_min_y)
          tile_geom = ogr.Geometry(ogr.wkbPolygon)
          tile_geom.AddGeometry(ring)

          # Check if tile intersects ROI
          if tile_geom.Intersects(roi_geom):
              subextents.append((subextent_min_x, subextent_max_x, subextent_min_y, subextent_max_y))

  return subextents


def get_IMG(wcs_url, fname, spatial_resolution=None):
  """

  Downloads and saves an image from a Web Coverage Service (WCS) URL.

  Parameters
  ----------

  wcs_url : str
    The URL of the Web Coverage Service (WCS) to download the image from.

  fname : str
    The output file path to save the downloaded image.

  spatial_resolution : float, optional
    The spatial resolution of the image in degrees.

  Returns
  -------
  str
    A string indicating the status of the download operation ("success" or "fail").
      
  """
  
  ds = gdal.Open(wcs_url)
  # Check if the dataset is successfully opened
  if ds is None:
      print(f"Failed to open WCS dataset from {wcs_url}")
      return "fail"
  else:
      output_format = "GTiff"
      # Create output dataset
      image_driver = gdal.GetDriverByName(output_format)

      if spatial_resolution:
          # Use specified spatial resolution if provided
          output_ds = gdal.Warp(
              fname,
              ds,
              xRes=spatial_resolution,
              yRes=spatial_resolution,
              resampleAlg="near",
              format=output_format
          )
      else:
          # Use default spatial resolution
          output_ds = image_driver.CreateCopy(fname, ds)

      # Close datasets
      ds = None
      output_ds = None

      print(f"Image saved to {fname}")
      return "success"


def build_mosaic(input_IMGs, output_mosaic):
  
  """
  Build a mosaic from a list of input images.

  Parameters
  ----------

  input_IMGs : list
    A list of input image paths to be mosaicked.

  output_mosaic : str
    The output mosaic image path.
      
  Returns
  -------
  None
  """

  if len(input_IMGs) == 1:
      # Only one tile, just copy it
      gdal.Translate(output_mosaic, input_IMGs[0])
      print(f"Mosaic created: {output_mosaic}")
      return

  # Warp options - set the format to GeoTIFF and compression
  warp_options = gdal.WarpOptions(
      format="GTiff",
      srcNodata=0,
      dstNodata=0,
      resampleAlg="near",
      creationOptions=["COMPRESS=LZW"]
  )

  # Mosaic the input files
  gdal.Warp(output_mosaic, input_IMGs, options=warp_options)
  print(f"Mosaic created: {output_mosaic}")


def get_url(extent, year):
  """
  Generate the WCS URL for the specified extent and year.

  Parameters
  ----------
  extent : tuple
    A tuple containing the minimum and maximum x and y coordinates of the extent.

  year : int
    The year of the NLCD data to be downloaded.

  Returns
  -------
  str
    The WCS URL for the specified extent and year.
  """
  
  min_x, max_x, min_y, max_y = extent
  coverage_id = f"NLCD_{year}_Land_Cover_L48"
  url = (
      f"https://www.mrlc.gov/geoserver/mrlc_display/{coverage_id}/wcs?"
      f"service=WCS&version=2.0.1&request=getcoverage&coverageid={coverage_id}"
      f"&subset=Lat({min_y},{max_y})&subset=Long({min_x},{max_x})"
      f"&SubsettingCRS=http://www.opengis.net/def/crs/EPSG/0/4326"
  )
  return url


def get_roi_geometry(args):
  """
  Get the region of interest (ROI) geometry based on the input arguments.

  Parameters
  ----------
  args : argparse.Namespace
    The parsed command-line arguments.

  Returns
  -------
  ogr.Geometry
    An OGR geometry representing the region of interest.
  """

  target_srs = osr.SpatialReference()
  target_srs.ImportFromEPSG(4326)
  target_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)

  if args.shapefile:
      # Open the shapefile
      shp_driver = ogr.GetDriverByName("ESRI Shapefile")
      dataset = shp_driver.Open(args.shapefile)

      if dataset is None:
          print(f"Failed to open shapefile: {args.shapefile}")
          sys.exit(1)

      layer = dataset.GetLayer()
      layer_srs = layer.GetSpatialRef()

      # Reproject to 4326 if necessary
      if not layer_srs.IsSame(target_srs):
          transform = osr.CoordinateTransformation(layer_srs, target_srs)

      # Merge all geometries into one
      roi_geom = ogr.Geometry(ogr.wkbMultiPolygon)
      for feature in layer:
          geom = feature.GetGeometryRef()
          if not layer_srs.IsSame(target_srs):
              geom.Transform(transform)
          roi_geom.AddGeometry(geom.Clone())
      roi_geom = roi_geom.UnionCascaded()
      dataset = None

  elif args.point:
      lat, lon = args.point
      point = ogr.Geometry(ogr.wkbPoint)
      point.AddPoint(lon, lat)
      # Buffer the point to create a small polygon (e.g., 0.01 degrees)
      buffer_distance = 0.0003  # degrees
      roi_geom = point.Buffer(buffer_distance)

  elif args.extent:
      minx, miny, maxx, maxy = args.extent
      ring = ogr.Geometry(ogr.wkbLinearRing)
      ring.AddPoint(minx, miny)
      ring.AddPoint(minx, maxy)
      ring.AddPoint(maxx, maxy)
      ring.AddPoint(maxx, miny)
      ring.AddPoint(minx, miny)
      roi_geom = ogr.Geometry(ogr.wkbPolygon)
      roi_geom.AddGeometry(ring)
  else:
      print("No valid ROI provided.")
      sys.exit(1)

  return roi_geom


if __name__ == "__main__":
  main()