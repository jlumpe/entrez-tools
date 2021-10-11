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


def one_or_none(items, msg: str ='Got more than one result'):
	"""Get only element of collection or None if empty, raise error if it contains more than one element."""
	if len(items) == 0:
		return None
	if len(items) == 1:
		return items[0]
	raise ValueError(msg)
