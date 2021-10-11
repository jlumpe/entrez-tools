

def get_esummary_result_json(data: dict) -> dict:
	"""Extract single result from ESummary JSON data."""
	result = data['result']
	uids = result['uids']
	assert len(uids) == 1, 'Data should contain exactly one result'
	return result[uids[0]]
