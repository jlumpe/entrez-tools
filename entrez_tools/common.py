import re


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
