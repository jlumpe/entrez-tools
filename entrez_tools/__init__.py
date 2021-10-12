"""
Tools to make working with the Entrez EUtils more bearable.
"""

__author__ = 'Jared Lumpe'
__email__ = 'mjlumpe@gmail.com'
__version__ = '0.1'


from .common import Database, DATABASES
from .elink import get_elink_map, get_elink_map_single
from .esearch import make_esearch_term, esearch_unique, esearch_accession
from .esummary import get_esummary_result_json
