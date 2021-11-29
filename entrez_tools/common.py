import re
from typing import List
from warnings import warn

from Bio.Entrez import Parser


def _check_xml_attrs(value, strict, path):
	if value.attributes:
		path_str = '/'.join(map(str, path))
		msg = f'{type(value)} at {path_str} has attributes.'
		if strict:
			raise ValueError(msg)
		else:
			warn(msg)

def xml_to_builtin(value, strict_attrs=False, _path=()):
	"""Recursively convert parsed XML value to builtin Python value.

	Result is similar to parsed JSON. Assumes XML elements do not have any attributes.

	Parameters
	----------
	value
		Value to convert.
	strict_attrs
		If True raise error when element has attributes, otherwise just warn.
	"""
	if isinstance(value, Parser.NoneElement):
		_check_xml_attrs(value, strict_attrs, _path)
		return None

	if isinstance(value, Parser.IntegerElement):
		_check_xml_attrs(value, strict_attrs, _path)
		return int(value)

	if isinstance(value, Parser.StringElement):
		_check_xml_attrs(value, strict_attrs, _path)
		return str(value)

	if isinstance(value, Parser.ListElement):
		_check_xml_attrs(value, strict_attrs, _path)
		return [
			xml_to_builtin(x, strict_attrs, (*_path, i))
			for i, x in enumerate(value)
		]

	if isinstance(value, Parser.DictionaryElement):
		_check_xml_attrs(value, strict_attrs, _path)
		return {
			k: xml_to_builtin(v, strict_attrs, (*_path, k))
			for k, v in value.items()
		}

	return value


def one_or_none(items, msg: str ='Got more than one result'):
	"""Get only element of collection or None if empty, raise error if it contains more than one element."""
	if len(items) == 0:
		return None
	if len(items) == 1:
		return items[0]
	raise ValueError(msg)


def lookslike_uid(s: str) -> bool:
	"""Check if a string looks like a UID (all digits)."""
	return re.fullmatch(r'\d+', s) is not None
