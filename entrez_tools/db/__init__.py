"""Tools for specific NCBI databases."""

from types import MappingProxyType
import typing as t
from attr import attrs, attrib


@attrs(frozen=True)
class Database:
	"""NCBI/Entrez database.

	Attributes
	----------
	name
		Database name.
	accession_search_field
		ESearch field name to search by accession number.
	accession_esummary_attr
		ESummary property name (JSON format) containing the accession number.
	"""
	name: str = attrib()
	accession_search_field: str = attrib()
	accession_esummary_attr: str = attrib()


DATABASES = MappingProxyType({
	'nuccore': Database(
		'nuccore',
		'Accession',
		'accessionversion',
	),
	'assembly': Database(
		'assembly',
		'Assembly Accession',
		'assemblyaccession',
	),
	# TODO...
})


DbArg = t.Union[str, Database]

def get_db(db: DbArg) -> Database:
	if isinstance(db, Database):
		return db
	if isinstance(db, str):
		return DATABASES[db]
	raise TypeError(f'Expected str or Database, got {db!r}')
