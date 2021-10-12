
def assembly_seq_url_from_ftppath(ftppath: str) -> str:
	"""Get FTP URL of assembly sequence file from FTP directory path."""
	return ftppath + '/' + ftppath.rsplit('/', 1)[1] + '_genomic.fna.gz'

def assembly_seq_url_from_esummary(esummary: dict, refseq: bool = True) -> str:
	"""Get FTP URL of assembly sequence file from ESummary JSON data."""
	attr = 'ftppath_refseq' if refseq else 'ftppath_genbank'
	return assembly_seq_url_from_ftppath(esummary[attr])
