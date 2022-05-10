from Bio.Blast import NCBIXML
from collections import defaultdict
class CovidVariant:
   def __init__(self, nm, aa, mt):
      self.name = nm
      self.protein = aa
      self.mutation = mt
      
def readFile(filename):
   f = open(filename)
   genome = filename.split(".",1)[0]
   curr_seq = []
   for line in f:
      line = line.strip()
      if len(line) != 0 and not line.startswith(">"): curr_seq.append(line)
   return ''.join(curr_seq)

def compareGenomes(variants, filename, query, sbjcr):
   # changes = defaultdict(lambda: defaultdict(list))
   mutation = defaultdict(lambda: 0)
   variantAA = defaultdict(lambda:0)
   
   f = open(filename,"r")
   item = next(NCBIXML.parse(f))
   previous = ["","","",""]
   for alignment in item.alignments:
      for hsp in alignment.hsps:
         if hsp.expect < 0.01:
            for num in range(len(hsp.query)):
               if hsp.query[num] != hsp.sbjct[num]:
                  mt = hsp.query[num] + str(num+1) + hsp.sbjct[num]
                  location = num
                  if previous[0] != "" and num == int(previous[1]):
                     location = mutation[previous[3]]
                     mt = previous[0] + hsp.query[num] + previous[1]+ previous[2]+ hsp.sbjct[num]
                     del mutation[previous[3]]
                  mutation[mt] = location
                  previous = [previous[0] + hsp.query[num], str(num+1), previous[2]+ hsp.sbjct[num], mt]
      variants[sbjcr].append(CovidVariant(query,rsAminoAcid, mutation))

#================================================#
genomes = {}
f = open("genomes.txt")
rs = readFile("references_sequences.fasta")
variants = defaultdict(list)

rsAminoAcid = defaultdict(lambda:0)
for aa in rs:
   rsAminoAcid[aa] +=1
reqSeg = CovidVariant("RegSeq",rsAminoAcid,None)


for line in f:
   elements = line.strip()
   if len(elements) != 0:
      elements = elements.split()
      filename, query, sbjcr = elements[0], elements[1], elements[2]
      compareGenomes(variants, filename, query, sbjcr)
   
for name, variants in variants.items():
   print(f"{name}'s S protein mutation compare to:")
   sameVariant = 2
   for i, v in enumerate(variants):
      otherName = v.name
      if v.name == name:
         otherName += " " + str(sameVariant)
         sameVariant+=1
      print(f"\t{i+1}. {otherName}: {[x for x in v.mutation.keys()]}")