from collections import namedtuple
# Imports for scraping miRBase
from bs4 import BeautifulSoup
import urllib2
import time
import re
# Imports for accessing Entrez via BioPython
import sys
from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW
Entrez.email = "rpeace@sce.carleton.ca"

# Extract all miRNA names and relative locations on the mirBase website
Mirna = namedtuple('Mirna', 'id chr start stop strand')
ggaPage = urllib2.urlopen("http://mirbase.org/cgi-bin/mirna_summary.pl?org=ath")
soup = BeautifulSoup(ggaPage)
# get IDs from page
ggaLinkSoup = soup.find_all(href=re.compile("\?acc="))[::2]
ids = [link.get_text() for link in ggaLinkSoup]
# Get chromosome, start, stop, and strand data from page
chrCols = soup.find_all(attrs={'class':'chrCol'})
chrs = [c.get_text() for c in chrCols[4::4]]
starts = [c.get_text() for c in chrCols[5::4]]
stops = [c.get_text() for c in chrCols[6::4]]
strands = [c.get_text() for c in chrCols[7::4]]
# Zip all the data so that we can process it more easily
ggaData = zip(ids, chrs, starts, stops, strands)

mirnaList = []
for datum in ggaData:
	if datum[1] != '':
		if True:
			mirnaList.append(Mirna(datum[0],
				datum[1],
				int(datum[2])-30,
				int(datum[3])+30,
				datum[4]))

print mirnaList
mirnaList = sorted(mirnaList, key=lambda m: m.chr)
chrList = []
for m in mirnaList:
	if m.chr not in chrList:
		chrList.append(m.chr)
if "X" in chrList:
	chrList.remove("X")

with open('data/ath.fasta', 'w') as outFile:
	print "CHRLIST: ", chrList
	for i in chrList:
		if i == 'chr5':
			searchTerm = "332002898"
		else:
			searchTerm = "Arabidopsis thaliana chromosome "+i[3]+" complete sequence PRJNA10719"
		# searchTerm = "Caenorhabditis elegans[organism] "+i+"[title] complete sequence bristol N2 chromosome "
		# for j in range(1,maxChr+1):
		# 	if j != i:
		# 		searchTerm += " NOT "+str(j)+"[title]"
		handle = Entrez.esearch(db="nucleotide", term=searchTerm)
		# handle = Entrez.esearch(db="nucleotide", term="(homo sapiens[Organism] chromosome "+str(i)+"[Title] grch37 primary reference assembly)NOT refseqgene NOT contig NOT patch NOT scaffold")
		record = Entrez.read(handle)
		gi = record["IdList"][0]
		print gi
		fetch_result = Entrez.efetch(db="nucleotide", id=gi, rettype="fasta", retmode="text")
		seqRecord = SeqIO.read(fetch_result, "fasta")
		for mirna in mirnaList:
			if mirna.chr == i:
				outFile.write(">"+mirna.id.replace('\n','')+'...'+str(mirna.chr)+':'+str(mirna.start)+'-'+str(mirna.stop)+'\n')
				if mirna.strand == '+':
					outFile.write(str(seqRecord.seq[mirna.start:mirna.stop].reverse_complement())+'\n')
				else:
					outFile.write(str(seqRecord.seq[mirna.start:mirna.stop])+'\n')