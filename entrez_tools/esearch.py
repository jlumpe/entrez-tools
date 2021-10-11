import typing as t

from Bio import Entrez
from Bio.Entrez.Parser import DictionaryElement

from .common import one_or_none, get_db, DbArg


def make_esearch_term(field: str, value) -> str:
	"""Make ESearch term string for field."""
	return f'({value}[{field}])'


def validate_esearch_result(result: DictionaryElement) -> None:
	"""Check ESearch result has no errors."""
	if 'ErrorList' in result:
		err = result['ErrorList']
		fieldnotfound = err.get('FieldNotFound', [])
		if fieldnotfound:
			raise ValueError(f'Invalid field(s): {fieldnotfound!r}')


def get_esearch_ids(result: DictionaryElement) -> t.List[str]:
	"""Get ID list from ESearch result."""
	return result['IdList']


def esearch_unique(db: str, field: str, value, **kw) -> t.Optional[str]:
	"""ESearch based on single field and get a unique result or None."""
	term = make_esearch_term(field, value)
	result = Entrez.read(Entrez.esearch(db, term=term, **kw))
	validate_esearch_result(result)
	ids = get_esearch_ids(result)
	return one_or_none(ids)


def esearch_accession(db: DbArg, acc: str, **kw) -> t.Optional[str]:
	"""ESearch query to look up UID by accession no."""
	field = get_db(db).accession_search_field
	return esearch_unique(db, field, acc, **kw)
