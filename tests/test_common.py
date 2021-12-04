import pytest

from entrez_tools.common import one_or_none, get_uid_list, lookslike_uid


def test_one_or_none():
	assert one_or_none([]) is None
	assert one_or_none(['foo']) == 'foo'
	with pytest.raises(ValueError):
		one_or_none(['foo', 'bar'])


def test_get_uid_list():
	uids = ['123', '456', '789']

	assert get_uid_list(uids) == uids
	assert get_uid_list([123, 456, 789]) == uids
	assert get_uid_list('123,456,789') == uids
	assert get_uid_list('123') == ['123']
	assert get_uid_list(123) == ['123']
	assert get_uid_list([]) == []


def test_lookslike_uid():
	assert lookslike_uid('123456')
	assert not lookslike_uid('GCF_000000.0')
	assert not lookslike_uid('')
