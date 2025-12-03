"""
S-102 Python Utilities

Tools for converting various datasets to S-102.
"""
from .api import S102File, S102Exception, PRODUCT_SPECIFICATION
from .utils import (
    create_s102,
    add_metadata,
    add_bathymetry_instance,
    add_data_from_arrays,
    update_metadata,
    write_data_file,
    from_arrays,
    from_arrays_with_metadata,
    from_gdal,
    from_bag,
    to_geotiff,
)

__all__ = ["api", "utils"]