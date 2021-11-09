"""Work with assembly database."""

import typing as t
import xml.etree.ElementTree as ET

import pandas as pd


def try_to_int(s: str):
	try:
		return int(s)
	except ValueError:
		return s


def seq_url_from_ftppath(ftppath: str) -> str:
	"""Get FTP URL of assembly sequence file from FTP directory path."""
	return ftppath + '/' + ftppath.rsplit('/', 1)[1] + '_genomic.fna.gz'

def seq_url_from_esummary(esummary: dict, refseq: bool = True) -> str:
	"""Get FTP URL of assembly sequence file from ESummary JSON data."""
	attr = 'ftppath_refseq' if refseq else 'ftppath_genbank'
	return seq_url_from_ftppath(esummary[attr])


def parse_summary_meta(s: str) -> ET.Element:
	"""Parse XML content of "meta" property in summary."""
	return ET.fromstring('<Meta>' + s + '</Meta>')

def format_summary_meta(meta: t.Union[str, ET.Element], as_int: bool = True) -> dict:
	"""Format parsed meta XML to dict."""
	if isinstance(meta, str):
		meta = parse_summary_meta(meta)

	out = dict()

	for child in meta:
		if child.tag == 'Stats':
			out.update(summary_meta_stats_dict(child, as_int=as_int))
		elif child.tag == 'FtpSites':
			pass  # TODO
		elif len(child) == 0 and not(child.attrib):
			out[child.tag] = child.text
		else:
			print(f'Don\'t know how to process element <{child.tag}>')

	return out

def summary_meta_stats_table(stats: ET.Element, as_int: bool = True) -> pd.DataFrame:
	"""Create table from assembly meta <Stats> element."""
	rows = [
		(stat.attrib['category'], stat.attrib['sequence_tag'], try_to_int(stat.text) if as_int else stat.text)
		for stat in stats
	]
	return pd.DataFrame.from_records(rows, columns=['category', 'sequence_tag', 'value'])

def summary_meta_stats_dict(stats: ET.Element, omit_all: bool = True, as_int: bool = True) -> dict:
	"""Create dictionary from assembly meta <Stats> element."""
	out = dict()

	for stat in stats:
		key = stat.attrib['category']
		tag = stat.attrib['sequence_tag']
		if not (omit_all and tag == 'all'):
			key += '_' + tag
		out[key] = try_to_int(stat.text) if as_int else stat.text

	return out
