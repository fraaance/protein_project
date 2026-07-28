################################################################
# The Reader script receives a FASTA file of .faa- or fna-format
# In case the RNA is not translated yet, the script will create a 
# translated protein sequence

# Finally, the script stores all the sequences of the received FASTA file
# as a protein (with its RNA seq. if available), storing all information of the FASTA
################################################################

from dataclasses import dataclass
import re

# FastaReader receives file path and returns protein name, organism, ID, sequence

class FastaReader: 
    def __init__(self, file_path):
        val, file_type = self.validate_path(file_path)
        self.prot_list = []

        if val:
            print("*" * 50)
            if file_type == "prot":
                print(" FASTA - Amino Acid Sequence".center(50, '*'))
                header, seq = self.read_file(file_path)
                accession, gene_name = re.split(r"\s+", header[0])
                gene_id = int(header[2].replace("GeneID=", "").rstrip("]"))
                accession=accession[1:]
                protein = Protein(accession=accession, 
                                gene_name=gene_name, 
                                organism=header[1].replace("organism=", ""), 
                                gene_id=gene_id, 
                                sequence=seq)
                print(f"Protein {gene_name} [GeneID={gene_id}]".ljust(50, " "))
                print(seq)
                self.prot_list.append(protein)
            
            if file_type == "dna":
                print(" FASTA - DNA Sequence ".center(50, '*'))
                print("*" * 50, "\n")
                print("Translating DNA into Protein".center(50, "*"))
                header, seq = self.read_file(file_path)
            
                #print(re.split(r"[:\-\s]+|(?<=c)(?=\d)", header[0]))
                accession, strand, start, end, gene_name = re.split(r"[:\-\s]+|(?<=c)(?=\d)", header[0])
                gene_id=int(header[2].replace("GeneID=", "").rstrip("]"))
                seq = decode_rna_to_protein(seq.replace("T", "U"))
                protein = Protein(accession="None",
                        gene_name=gene_name, 
                        organism=header[1].replace("organism=", "").rstrip("]"), 
                        gene_id=gene_id, 
                        sequence=seq)
                print(f"Protein {gene_name} [GeneID={gene_id}]".ljust(50, " "))
                print(seq)
                #print(f"\nsequence {gene_name} [GeneID={gene_id}] successfully read".center(50, ' '))
                self.prot_list.append(protein)
        

    def validate_path(self, file_path):
        if (file_path.endswith(".faa")):
            return True, "prot"
        elif (file_path.endswith(".fna")):
            return True, "dna"
        return False

    def read_file(self, file_path):
        with open(file_path) as f:
            lines = f.read().splitlines()
            header = re.split(r"\s\[+", lines[0])
            seq = "".join(lines[1:])
        return header, seq

    def print_out(self, input):
        print("*" * 50)
        print(input.accession, input.gene_name, input.gene_id, input.organism)
        print(len(input.sequence), input.sequence)

def decode_rna_to_protein(sequence):
    # https://github.com/T101J/Translating_RNA_to_Protein.git
    # source 

    rna_codons = {
        "UUU" : "F", "CUU" : "L", "AUU" : "I", "GUU" : "V",
        "UUC" : "F", "CUC" : "L", "AUC" : "I", "GUC" : "V",
        "UUA" : "L", "CUA" : "L", "AUA" : "I", "GUA" : "V",
        "UUG" : "L", "CUG" : "L", "AUG" : "M", "GUG" : "V",
        "UCU" : "S", "CCU" : "P", "ACU" : "T", "GCU" : "A",
        "UCC" : "S", "CCC" : "P", "ACC" : "T", "GCC" : "A",
        "UCA" : "S", "CCA" : "P", "ACA" : "T", "GCA" : "A",
        "UCG" : "S", "CCG" : "P", "ACG" : "T", "GCG" : "A",
        "UAU" : "Y", "CAU" : "H", "AAU" : "N", "GAU" : "D",
        "UAC" : "Y", "CAC" : "H", "AAC" : "N", "GAC" : "D",
        "UAA" : "STOP", "CAA" : "Q", "AAA" : "K", "GAA" : "E",
        "UAG" : "STOP", "CAG" : "Q", "AAG" : "K", "GAG" : "E",
        "UGU" : "C", "CGU" : "R", "AGU" : "S", "GGU" : "G",
        "UGC" : "C", "CGC" : "R", "AGC" : "S", "GGC" : "G",
        "UGA" : "STOP", "CGA" : "R", "AGA" : "R", "GGA" : "G",
        "UGG" : "W", "CGG" : "R", "AGG" : "R", "GGG" : "G" 
        }
    
    protein_seq = ""
    #print(len(sequence)/3)
    for i in range(0, len(sequence) - (3 + len(sequence)%3), 3):
        codon = rna_codons[sequence[i:i+3]]
        if codon == "STOP":
            break
        protein_seq += codon
    return protein_seq
        
@dataclass
class RNA:
    accession: str
    #start: int
    #end: int 
    #strand: str
    gene_name: str 
    organism: str
    gene_id: int 
    #chromosome: str | None
    sequence: str 
        
@dataclass
class Protein: 
    accession: str 
    gene_name: str 
    organism: str
    gene_id: int
    sequence: str 




gene_file = FastaReader("/Users/franzweisel/Downloads/project/nadE_NAD_synthetase/data/gene.fna")
protein_file = FastaReader("/Users/franzweisel/Downloads/project/nadE_NAD_synthetase/data/protein.faa")

#test_file = FastaReader("/Users/franzweisel/Documents/Systembio/Cutibacterium_granulosum_TM11/cuti_TM11_translated_cds.faa")
test2_file = FastaReader("/Users/franzweisel/Downloads/ncbi_dataset-5/ncbi_dataset/data/protein.faa")