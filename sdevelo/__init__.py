import os
from time import gmtime, strftime

from ._sim import SimData
from ._config import Config
from ._model import SDENN
from ._pl import plot_streamline, plot_latent_time, plot_noise_histogram, plot_gene_scatter, plot_subset
from ._infer import infer_gene_correlations

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

try:
    from setuptools_scm import get_version
    __version__ = get_version(root="..", relative_to=__file__)
    del get_version
    
except (LookupError, ImportError):
    try:
        from importlib_metadata import version 
    except:
        from importlib.metadata import version 
    __version__ = version(__name__)
    del version

print (f'(Working on SDEvelo {__version__})')
print (strftime("%Y-%m-%d %H:%M:%S", gmtime()))


