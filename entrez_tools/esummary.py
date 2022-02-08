import json

from Bio import Entrez

from .db import get_db, DbArg
from .common import get_uid_list


def get_esummary_json(data: dict, uids=None) -> dict:
	"""Validate parsed ESummary JSON response and get values.

	Parameters
	----------
	data
		Parsed JSON response.
	uids
		Assert that result contains these UIDs.

	Returns
	-------
	dict
		Dict mapping UIDs to summaries.
	"""
	result = data['result']

	if uids is not None:
		uids = get_uid_list(uids)
		assert set(uids) == set(result['uids'])

	out = dict()

	for uid in result['uids']:
		summary = result[uid]
		assert not summary.get('error')
		out[uid] = summary

	return out


def esummary_json(db: str, ids, strict=True):
	"""Fetch ESummary data in JSON format.

	Parameters
	----------
	db
	ids
	strict
		Validate result.

	Returns
	-------
	dict
		Dict mapping UIDs to summaries.
	"""
	# Current version of BioPython doesn't accept a list
	ids_str = ','.join(get_uid_list(ids))
	with Entrez.esummary(db=db, id=ids_str, retmode='json') as response:
		data = json.load(response)
	return get_esummary_json(data, uids=ids if strict else None)


def get_esumary_accession(db: DbArg, summary: dict) -> str:
	"""Get accession no. from ESummary JSON dict."""
	return summary[get_db(db).accession_esummary_attr]
