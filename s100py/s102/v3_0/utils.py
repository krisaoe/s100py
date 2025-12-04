""" Functions to create S102 data from other sources

If the utils module is run as __main__ then it will run  :func:`from_bag` or :func:`from_gdal` on the filename given in the command line arguments.

"""
# this also works but doesn't put parenthesis in the html   :any:`from_bag`
# fully qualified would also work   :func:`s100py.s102.utils.from_gdal`

import os
import sys

if getattr(sys, 'frozen', False):
    # in a frozen exe os.pyc is at the root level
    proj_db_path = os.path.join(os.path.dirname(os.__file__), "library\\share\\proj")
    # print("frozen exe", proj_db_path)
    os.environ["PROJ_LIB"] = proj_db_path
    # print(os.listdir(proj_db_path))  # os.listdir(os.path.dirname(os.__file__)))
else:
    pass
    # print("running as script")

import argparse
import datetime
import logging
import warnings
from typing import Any, Dict, Mapping, Optional, Union
from xml.etree import ElementTree as et
import tkinter as tk
from tkinter import filedialog, messagebox
import numpy
from osgeo import gdal, osr
import h5py
import re

try:
    from matplotlib import pyplot
    from matplotlib.colors import LinearSegmentedColormap, ListedColormap, BoundaryNorm
    from matplotlib.colorbar import ColorbarBase
except:
    if not getattr(sys, 'frozen', False):  # we expect the frozen exe to not have matplotlib
        print("matplotlib.pyplot failed to import, plotting will not work")

from s100py.s102.v3_0.api import DEPTH, UNCERTAINTY, S102File, S102Exception
from ...s1xx import s1xx_sequence


Number = Union[int, float]
ArrayLike = Union[numpy.ndarray, s1xx_sequence]
MetadataMapping = Mapping[str, Any]
GridProps = Mapping[str, Number]


def create_s102(output_path: str, data_coding_format: int = 2) -> S102File:
    """
    Create S-102 v3.0.0 file.

    Parameters
    ----------
    output_path : str
        Path to the HDF5 file to create.
    data_coding_format : int
        Data coding format (DCF). Must be 2 (regular grid) for S-102.

    Returns
    -------
    S102File
        Open S-102 file object.
    """
    if data_coding_format != 2:
        raise S102Exception(f"S102 only supports data_coding_format=2 (regular grid), got {data_coding_format}")

    return S102File.create_s102(output_path)


def add_metadata(
    metadata: MetadataMapping,
    data_file: S102File,
) -> S102File:
    """
    Populate root / Group_F metadata for an S-102 v3.0.0 dataset.

    `metadata` is a pythonic dict of values; this function maps them to
    S102Root + Group_F attributes.

    Typical metadata keys you might support (align with S-102 3.0.0):

      - horizontalCRS  (EPSG integer)
      - verticalCS
      - verticalDatumReference
      - verticalCoordinateBase
      - verticalDatum
      - epoch
      - verticalDatumEpoch
      - dataDynamicity
      - issueDateTime
      - geographicIdentifier
      - datasetID, editionNumber, updateNumber
      - depthUncertainty, horizontalPositionUncertainty, verticalUncertainty
      - datasetDeliveryInterval, timeUncertainty
      - producingAgency, etc.

    This is analogous to s104.utils.add_metadata().
    """
    root = data_file.root

    # --- horizontal / vertical reference frames ---
    if "horizontalCRS" in metadata and hasattr(root, "horizontal_crs"):
        root.horizontal_crs = int(metadata["horizontalCRS"])

    if "verticalCS" in metadata and hasattr(root, "vertical_cs"):
        root.vertical_cs = int(metadata["verticalCS"])

    if "verticalDatumReference" in metadata and hasattr(root, "vertical_datum_reference"):
        root.vertical_datum_reference = int(metadata["verticalDatumReference"])

    if "verticalCoordinateBase" in metadata and hasattr(root, "vertical_coordinate_base"):
        root.vertical_coordinate_base = int(metadata["verticalCoordinateBase"])

    if "verticalDatum" in metadata and hasattr(root, "vertical_datum"):
        root.vertical_datum = int(metadata["verticalDatum"])

    if "epoch" in metadata and hasattr(root, "epoch"):
        root.epoch = str(metadata["epoch"])

    if "verticalDatumEpoch" in metadata and hasattr(root, "vertical_datum_epoch"):
        root.vertical_datum_epoch = str(metadata["verticalDatumEpoch"])

    # --- identification / edition ---
    if "datasetID" in metadata and hasattr(root, "dataset_id"):
        root.dataset_id = str(metadata["datasetID"])

    if "editionNumber" in metadata and hasattr(root, "edition_number"):
        root.edition_number = int(metadata["editionNumber"])

    if "updateNumber" in metadata and hasattr(root, "update_number"):
        root.update_number = int(metadata["updateNumber"])

    if "geographicIdentifier" in metadata and hasattr(root, "geographic_identifier"):
        root.geographic_identifier = str(metadata["geographicIdentifier"])

    # --- time / dynamicity ---
    if "dataDynamicity" in metadata and hasattr(root, "data_dynamicity"):
        root.data_dynamicity = int(metadata["dataDynamicity"])

    if "issueDateTime" in metadata and hasattr(root, "issue_date_time"):
        val = metadata["issueDateTime"]
        if isinstance(val, datetime.datetime):
            root.issue_date_time = val
        else:
            # accept ISO 8601 string, etc.
            root.issue_date_time = datetime.datetime.fromisoformat(str(val))

    if "datasetDeliveryInterval" in metadata and hasattr(root, "dataset_delivery_interval"):
        root.dataset_delivery_interval = str(metadata["datasetDeliveryInterval"])

    if "timeUncertainty" in metadata and hasattr(root, "time_uncertainty"):
        root.time_uncertainty = float(metadata["timeUncertainty"])

    # --- uncertainties (often Group_F / Feature Attribute Table) ---
    group_f = getattr(root, "group_f", None)
    if group_f is not None:
        # This assumes your Group_F has feature attribute records where
        # you can set default uncertainty values. How you do that will
        # depend on the v3.0.0 FC and the generated classes.
        if "depthUncertainty" in metadata and hasattr(group_f, "default_depth_uncertainty"):
            group_f.default_depth_uncertainty = float(metadata["depthUncertainty"])
        if "horizontalPositionUncertainty" in metadata and hasattr(
            group_f, "default_horizontal_position_uncertainty"
        ):
            group_f.default_horizontal_position_uncertainty = float(
                metadata["horizontalPositionUncertainty"]
            )
        if "verticalUncertainty" in metadata and hasattr(group_f, "default_vertical_uncertainty"):
            group_f.default_vertical_uncertainty = float(metadata["verticalUncertainty"])

    # You can also wire in producingAgency etc. here when you know the exact names.
    return data_file



def add_bathymetry_instance(
    data_file: S102File,
    data_coding_format: int = 2,
) -> S102File:
    """
    Ensure that a BathymetryCoverage instance (and its group) exist.

    This is the S-102 analogue of s104.add_water_level_instance().

    - Creates root.bathymetry_coverage if needed.
    - Ensures there is at least one BathymetryCoverage feature instance.
    - Ensures that instance has at least one BathymetryGroup for DCF2.
    """
    if data_coding_format != 2:
        raise S102Exception(f"S102 only supports data_coding_format=2 (regular grid), got {data_coding_format}")

    root = data_file.root

    container = getattr(root, "bathymetry_coverage", None)
    if container is None:
        raise S102Exception(
            "root.bathymetry_coverage is not initialised; "
            "check S102Root API or create via create_s102()."
        )

    # --- Instance list ---
    instances = getattr(container, "bathymetry_coverage", None)
    if instances is None:
        raise S102Exception(
            "container.bathymetry_coverage list missing; "
            "check S102 API generation."
        )

    if len(instances) == 0:
        instance = instances.append_new_item()
    else:
        instance = instances[0]

    groups = getattr(instance, "bathymetry_group", None)
    if groups is None:
        raise S102Exception(
            "instance.bathymetry_group list missing; "
            "check S102 API generation."
        )

    if len(groups) == 0:
        groups.append_new_item()

    return data_file


def add_data_from_arrays(
    depth: ArrayLike,
    grid_properties: Dict[str, Number],
    data_file: S102File,
    uncertainty: Optional[ArrayLike] = None,
    nodata_value: Optional[Number] = None,
    z_positive_down: bool = True,
) -> S102File:
    """
    Insert depth (and optionally uncertainty) arrays into an S-102 dataset (DCF2),
    and update the BathymetryCoverage instance geometry.

    Parameters
    ----------
    depth : ArrayLike
        2D numpy array with shape [ny, nx].
    grid_properties : dict
        {
          "minx": westBound (x),
          "maxx": eastBound (x),
          "miny": southBound (y),
          "maxy": northBound (y),
          "cellsize_x": dx,
          "cellsize_y": dy,
          "nx": number of columns,
          "ny": number of rows,
        }
    data_file : S102File
        The S102 file to populate.
    uncertainty : ArrayLike, optional
        2D numpy array with same shape as depth. If None, fills with nodata.
    nodata_value : Number, optional
        Value to use for missing data. Defaults to 1000000.0.
    z_positive_down : bool
        If False, depth values are negated.
    """
    ny, nx = depth.shape
    if int(grid_properties["nx"]) != nx or int(grid_properties["ny"]) != ny:
        raise S102Exception(
            "Grid properties (nx, ny) do not match array shape: "
            f"({grid_properties['nx']}, {grid_properties['ny']}) vs {depth.shape}"
        )

    if uncertainty is not None and depth.shape != uncertainty.shape:
        raise S102Exception(
            f"Depth and uncertainty arrays must have the same shape, "
            f"got {depth.shape!r} and {uncertainty.shape!r}"
        )

    fill_value = numpy.float32(nodata_value if nodata_value is not None else 1000000.0)

    depth_arr = numpy.array(depth, dtype="float32", copy=True)
    if not z_positive_down:
        depth_arr *= -1.0
    depth_arr = numpy.where(numpy.isnan(depth_arr), fill_value, depth_arr)

    if uncertainty is not None:
        uncert_arr = numpy.array(uncertainty, dtype="float32", copy=True)
        uncert_arr = numpy.where(numpy.isnan(uncert_arr), fill_value, uncert_arr)
    else:
        uncert_arr = numpy.full(depth.shape, fill_value, dtype="float32")

    minx = float(grid_properties["minx"])
    maxx = float(grid_properties["maxx"])
    miny = float(grid_properties["miny"])
    maxy = float(grid_properties["maxy"])
    dx = float(grid_properties["cellsize_x"])
    dy = float(grid_properties["cellsize_y"])

    root = data_file.root

    # --- Navigate to instance + group ---
    container = root.bathymetry_coverage
    instances = container.bathymetry_coverage
    if not instances:
        raise S102Exception(
            "No BathymetryCoverage instances present; "
            "call add_bathymetry_instance() first."
        )
    instance = instances[0]

    groups = instance.bathymetry_group
    if not groups:
        raise S102Exception(
            "No BathymetryGroup present; call add_bathymetry_instance() first."
        )
    group = groups[0]

    # --- Instance-level geometry (DCF2) ---
    # Names here are based on the public S-102 examples; adjust to your API.
    instance.grid_origin_longitude = minx
    instance.grid_origin_latitude = miny
    instance.grid_spacing_longitudinal = dx
    instance.grid_spacing_latitudinal = dy

    instance.num_points_longitudinal = nx
    instance.num_points_latitudinal = ny

    instance.west_bound_longitude = minx
    instance.east_bound_longitude = maxx
    instance.south_bound_latitude = miny
    instance.north_bound_latitude = maxy

    # If the group has its own extents/size, mirror them there too:
    if hasattr(group, "num_points_longitudinal"):
        group.num_points_longitudinal = nx
    if hasattr(group, "num_points_latitudinal"):
        group.num_points_latitudinal = ny

    # --- Write arrays via S102File convenience accessors ---
    data_file.depth[:] = depth_arr
    data_file.uncertainty[:] = uncert_arr

    if nodata_value is not None:
        if hasattr(data_file.depth, "no_data_value"):
            data_file.depth.no_data_value = numpy.float32(nodata_value)
        if hasattr(data_file.uncertainty, "no_data_value"):
            data_file.uncertainty.no_data_value = numpy.float32(nodata_value)

    return data_file


def update_metadata(
    data_file: S102File,
    grid_properties: GridProps,
    metadata: Optional[MetadataMapping] = None,
    depth: Optional[ArrayLike] = None,
) -> S102File:
    """
    Update root-level metadata that depends on the actual coverage:

      - west/east/south/north bounds (from grid_properties)
      - minimumDepth/maximumDepth (from depth array if provided, or from file)
      - optional overrides from metadata dict

    This is analogous to s104.utils.update_metadata().
    """
    root = data_file.root

    minx = float(grid_properties["minx"])
    maxx = float(grid_properties["maxx"])
    miny = float(grid_properties["miny"])
    maxy = float(grid_properties["maxy"])

    # horizontal extents at root
    if hasattr(root, "west_bound_longitude"):
        root.west_bound_longitude = minx
    if hasattr(root, "east_bound_longitude"):
        root.east_bound_longitude = maxx
    if hasattr(root, "south_bound_latitude"):
        root.south_bound_latitude = miny
    if hasattr(root, "north_bound_latitude"):
        root.north_bound_latitude = maxy

    # vertical extents from depth values
    if depth is None:
        try:
            depth = numpy.array(data_file.depth[:], dtype="float32")
        except Exception:
            depth = None

    if depth is not None:
        depth_arr = numpy.array(depth, dtype="float32", copy=False)
        finite = numpy.isfinite(depth_arr)
        if numpy.any(finite):
            min_z = float(depth_arr[finite].min())
            max_z = float(depth_arr[finite].max())

            if hasattr(root, "minimum_depth"):
                root.minimum_depth = min_z
            if hasattr(root, "maximum_depth"):
                root.maximum_depth = max_z

    # Optional explicit overrides from metadata
    if metadata is not None:
        if "minimumDepth" in metadata and hasattr(root, "minimum_depth"):
            root.minimum_depth = float(metadata["minimumDepth"])
        if "maximumDepth" in metadata and hasattr(root, "maximum_depth"):
            root.maximum_depth = float(metadata["maximumDepth"])

        if "westBoundLongitude" in metadata and hasattr(root, "west_bound_longitude"):
            root.west_bound_longitude = float(metadata["westBoundLongitude"])
        if "eastBoundLongitude" in metadata and hasattr(root, "east_bound_longitude"):
            root.east_bound_longitude = float(metadata["eastBoundLongitude"])
        if "southBoundLatitude" in metadata and hasattr(root, "south_bound_latitude"):
            root.south_bound_latitude = float(metadata["southBoundLatitude"])
        if "northBoundLatitude" in metadata and hasattr(root, "north_bound_latitude"):
            root.north_bound_latitude = float(metadata["northBoundLatitude"])

    return data_file


def write_data_file(data_file: S102File) -> None:
    """
    Write and close an S-102 v3.0.0 file.

    Mirrors s104.utils.write_data_file().
    """
    data_file.write()
    data_file.close()



__all__ = ['plot_depth_using_h5py', 'create_s102', 'from_arrays', 'from_arrays_with_metadata',
           'from_gdal', 'from_bag', 'get_valid_epsg', 'to_geotiff']

gco = "{http://www.isotc211.org/2005/gco}"

# @todo create a friendly name mapping to s102 nested location, then add s102 functions for "to dictionary" and "from dictionary" to api
#   that would make these functions easily invertable

r"""
from s100py.s102 import utils
navo_name = r"G:\Data\S102 Data\LA_LB_Area_GEO_reprojected.bag.navo_%d.h5"
utils.plot_depth_using_h5py(navo_name)
bag_names = [r"G:\Data\S102 Data\NBS_US5NYCAH_20200430.bag", r"G:\Data\S102 Data\NBS_US5NYCBH_20200429.bag", r"G:\Data\S102 Data\LA_LB_Area_GEO_reprojected.bag"]
fout = utils.from_bag(bag_name, bag_name+r".noaa.h5")
for bag_name in bag_names:
    fout = utils.from_bag(bag_name, bag_name+r".noaa.h5")
    fout.close()
    utils.plot_depth_using_h5py(bag_name + r".noaa.h5")
"""

def plot_depth_using_h5py(filename, enc_color=False):
    # filename = r"G:\Data\S102 Data\GlenS102Test\102USA15NYCAH200430.H5"
    # h5py.File(r"G:\Data\S102 Data\LA_LB_Area_GEO_reprojected.bag_%d.h5", mode="r", driver="family", memb_size=681574400)
    try:
        f = h5py.File(filename, mode="r", driver="family")
    except OSError as e:
        # see if the member size isn't right and then fall back to standard opening
        # OSError: Unable to open file (Family member size should be 681574400.  But the size from file access property is 2147483647)
        try:
            error_string = str(e)
            m = re.search(r'Family member size should be\s+(\d+)', error_string)
            if m:
                sz = int(m[1])
            f = h5py.File(filename, mode="r", driver="family", memb_size=int(sz))
        except:
            f = h5py.File(filename, mode="r")
    fill_val = 1000000
    try:
        d = f["BathymetryCoverage/BathymetryCoverage.01/Group.001/values"]['depth']
    except KeyError:
        try:
            d = f["BathymetryCoverage/BathymetryCoverage.001/Group.001/values"]['depth']
        except KeyError:
            d = f["SurfaceCurrent/SurfaceCurrent.01/Group_001/values"]['surfaceCurrentSpeed']
            fill_val = -9999
    d[d==fill_val] = numpy.nan

    # ud = numpy.flipud(d)
    if enc_color:
        colors = numpy.array([(255, 255, 255), (201, 237, 252), (167, 217, 251), (130, 202, 255), (97, 180, 255)]) / 255
        # bins = [-1000, -9.1, -5.4, -3.6, -1.8, ]
        bins = [-1000, -30, -22, -15, -12]
        norm = BoundaryNorm(bins, colors.shape[0])
        cmap = ListedColormap(colors)
        pyplot.imshow(d, interpolation='nearest', cmap=cmap, norm=norm)
    else:
        pyplot.imshow(d)
    pyplot.gca().invert_yaxis()
    pyplot.show()
    # cm = LinearSegmentedColormap.from_list("enc", colors, N=4)
    # im = pyplot.imshow(ud, cmap=cm, interpolation='nearest')
    # pyplot.colorbar(im)


def browse_files(question):
    # using tkinter since it is built in to python and smaller to distribute than PySide2 or wxPython in an executable
    root = tk.Tk()
    root.withdraw()
    # root.filename = tkFileDialog.askopenfilename(initialdir="/", title="Select file",
    #                                              filetypes=(("jpeg files", "*.jpg"), ("all files", "*.*")))
    file_path = filedialog.askopenfilename(title=question)
    return file_path


def bool_question(question, title="", icon="warning"):
    root = tk.Tk()
    root.withdraw()
    result = messagebox.askquestion(title, question, icon=icon)
    return result == 'yes'


def make_parser():
    parser = argparse.ArgumentParser(description='Convert a georeferenced file to S102')
    parser.add_argument("-?", "--show_help", action="store_true", help="show this help message and exit")
    parser.add_argument("-i", "--input_filename", help="full path to the file to be processed")
    parser.add_argument("-o", "--output_filename", help="output filename, default is same name as input with .h5 appended")
    parser.add_argument("-r", "--res", help="Resolution.  If the input file is a BAG then use attempt to use the given resolution" )
    return parser


def to_geotiff(input_path, output_path):
    s102_data = S102File(input_path)
    s102_data.to_geotiff(output_path)


from_arrays = S102File.from_arrays
from_arrays_with_metadata = S102File.from_arrays_with_metadata
from_gdal = S102File.from_raster
from_raster = S102File.from_raster
from_bag = S102File.from_bag
get_valid_epsg = S102File.get_valid_epsg

if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args()
    if args.show_help:
        parser.print_help()
        sys.exit()

    if not args.input_filename:
        path = browse_files('Choose file to convert to S102')
        if path:
            args.input_filename = path

    if args.input_filename:
        # get the output file path from the command line or the user if it wasn't specified
        if args.output_filename:
            output_name = args.output_filename
        else:
            output_name = args.input_filename + ".h5"
            if os.path.exists(output_name):
                if bool_question(f"{output_name} already exists, overwrite?", "Overwrite File"):
                    os.remove(output_name)

        # check if the data is a bag and should be sent to the from_bag function or just raster and send to from_gdal
        ds = gdal.Open(args.input_filename)
        drv = ds.GetDriver()
        if drv.GetDescription() == "BAG":
            from_bag(ds, output_name)
        else:
            from_gdal(ds, output_name)
