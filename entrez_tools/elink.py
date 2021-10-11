import typing as t

from Bio.Entrez.Parser import DictionaryElement

from .common import one_or_none


def validate_elink_item(item: DictionaryElement) -> None:
	assert not item['ERROR']


def get_elink_item_to_ids(item: DictionaryElement) -> t.List[str]:
	linksetdb = item['LinkSetDb']

	if len(linksetdb) == 0:
		return []

	assert len(linksetdb) == 1
	linkset = linksetdb[0]

	return [x['Id'] for x in linkset['Link']]


def get_elink_item_pair(item: DictionaryElement) -> t.Tuple[str, t.List[str]]:
	"""Get a (from_id, to_ids) pair from an ELink result item."""
	validate_elink_item(item)

	assert len(item['IdList']) == 1
	from_id = item['IdList'][0]

	to_ids = get_elink_item_to_ids(item)

	return from_id, to_ids


def get_elink_item_pair_single(item: DictionaryElement) -> t.Tuple[str, str]:
	"""Get a (from_id, to_id) pair from an ELink result item."""
	from_id, to_ids = get_elink_item_pair(item)
	to_id = one_or_none(to_ids, f'Got more than link from ID {from_id}')
	return from_id, to_id


def get_elink_map(items: t.Iterable[DictionaryElement]) -> t.Dict[str, t.List[str]]:
	"""Get dictionary mapping from_id -> to_ids for ELink result."""
	return dict(get_elink_item_pair(item) for item in items)


def get_elink_map_single(items: t.Iterable[DictionaryElement]) -> t.Dict[str, str]:
	"""Get dictionary mapping from_id -> to_id for ELink result."""
	return dict(get_elink_item_pair_single(item) for item in items)
