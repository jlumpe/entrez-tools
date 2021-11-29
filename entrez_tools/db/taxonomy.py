from typing import Dict, Optional, Iterable

from Bio import Entrez

from entrez_tools.common import xml_to_builtin


def fetch_taxonomy_tree(taxids: Iterable[str],
                        taxa: Optional[Dict[str, Dict]] = None,
                        count: int = 200,
                        ) -> Dict[str, Dict]:
	"""Download complete taxonomy tree.

	Parameters
	----------
	taxids
		UIDs of leaf taxa to be included.
	taxa
		Dictionary mapping UIDs to taxonomy data already downloaded.
	count
		Number of results to fetch per request.

	Returns
	-------
	Dict[str, Dict]
		Dict mapping UIDs to ESummary data in JSON format.
	"""
	if taxa is None:
		taxa = dict()

	to_download = set(taxids) - taxa.keys()

	for taxid, taxon in taxa.items():
		parent = taxon['ParentTaxId']
		if parent != '0' and parent not in taxa:
			to_download.add(parent)

	while to_download:
		ids_str = ','.join(list(to_download)[:count])
		with Entrez.efetch(db='taxonomy', id=ids_str) as f:
			result = Entrez.read(f)

		for taxon in result:
			taxon = xml_to_builtin(taxon)
			taxid = taxon['TaxId']
			parent = taxon['ParentTaxId']

			assert taxid not in taxa
			taxa[taxid] = taxon
			to_download.remove(taxid)

			if parent != '0' and parent not in taxa:
				to_download.add(parent)

	return taxa
