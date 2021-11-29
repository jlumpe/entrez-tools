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
	has_accession
		Whether this database's entries have an accession # in addition to a UID.
	accession_search_field
		ESearch field name to search by accession number.
	accession_esummary_attr
		ESummary property name (JSON format) containing the accession number.
	"""
	name: str = attrib()
	has_accession: bool = attrib()
	accession_search_field: t.Optional[str] = attrib()
	accession_esummary_attr: t.Optional[str] = attrib()


DATABASES = MappingProxyType({
	'nuccore': Database(
		'nuccore',
		True,
		'Accession',
		'accessionversion',
	),
	'assembly': Database(
		'assembly',
		True,
		'Assembly Accession',
		'assemblyaccession',
	),
	'taxonomy': Database(
		'assembly',
		False,
		None,
		None,
	),
})


DbArg = t.Union[str, Database]

def get_db(db: DbArg) -> Database:
	if isinstance(db, Database):
		return db
	if isinstance(db, str):
		return DATABASES[db]
	raise TypeError(f'Expected str or Database, got {db!r}')
