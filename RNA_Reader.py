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
            content_list = self.read_file(file_path)
            for content in content_list:
                print("*" * 50)
                if file_type == "prot":
                    print(" FASTA - Amino Acid Sequence recognized ".center(50, '*'))
                    header, seq = content
                    #accession, gene_name = re.split(r"\s+", header[0])
                    #gene_id = int(header[2].replace("GeneID=", "").rstrip("]"))
                    #accession=accession[1:]
                    accession = None
                    gene_name = None
                    organism = None
                    gene_id = None

                    protein = Protein(header=header,
                                    accession=accession, 
                                    gene_name=gene_name, 
                                    organism=organism,#header[1].replace("organism=", ""), 
                                    gene_id=gene_id, 
                                    sequence=seq)
                    print(f"Protein {gene_name} [GeneID={gene_id}]".ljust(50, " "))
                    print(seq)
                    self.prot_list.append(protein)
                
                if file_type == "dna":
                    print(" FASTA - DNA Sequence recognized".center(50, '*'))
                    print("Translating DNA into Protein".center(50, "*"))
                    header, seq = content
                
                    #print(re.split(r"[:\-\s]+|(?<=c)(?=\d)", header[0]))
                    #accession, strand, start, end, gene_name = re.split(r"[:\-\s]+|(?<=c)(?=\d)", header[0])
                    #gene_id=int(header[2].replace("GeneID=", "").rstrip("]"))
                    seq = decode_rna_to_protein(seq.replace("T", "U"))
                    accession = None
                    gene_name = None
                    organism = None
                    gene_id = None

                    protein = Protein(header=header,
                            accession="None",
                            gene_name=gene_name, 
                            organism=organism, #header[1].replace("organism=", "").rstrip("]"), 
                            gene_id=gene_id, 
                            sequence=seq)
                    print(f"Protein {gene_name} [GeneID={gene_id}]".ljust(50, " "))
                    print(seq)
                    #print(f"\nsequence {gene_name} [GeneID={gene_id}] successfully read".center(50, ' '))
                    self.prot_list.append(protein)
                    
        print((f" Number of found Sequences: {len(self.prot_list)} ".center(50, "*")))

    def validate_path(self, file_path):
        if (file_path.endswith(".faa")):
            return True, "prot"
        elif (file_path.endswith(".fna")):
            return True, "dna"
        return False

    def read_file(self, file_path):
        with open(file_path) as f:
            liste = re.split(">", f.read())[1:]
            content_list = []
            print(liste[0])

            for e in liste:
                lines = e.splitlines()
                header = lines[0]
                seq = "".join(lines[1:])
                print("header: ", header)

                content_list.append((header, seq))
                #print(re.split(header))
                
        #with open(file_path) as f:
        #    lines = f.read().splitlines()
        #    header = re.split(r"\s\[+", lines[0])
        #    seq = "".join(lines[1:])
        return content_list

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
    header: str
    accession: str | None = None
    gene_name: str | None = None
    organism: str | None = None
    gene_id: int | None = None
    sequence: str = ""





####################################################
# Program calls
####################################################
gene_file = FastaReader("/Users/franzweisel/Downloads/project/nadE_NAD_synthetase/data/gene.fna")

protein_file = FastaReader("/Users/franzweisel/Downloads/project/nadE_NAD_synthetase/data/protein.faa")
test_file = FastaReader("/Users/franzweisel/Documents/Systembio/Cutibacterium_granulosum_TM11/cuti_TM11_translated_cds.faa")
test2_file = FastaReader("/Users/franzweisel/Downloads/ncbi_dataset-5/ncbi_dataset/data/protein.faa")

test3_file = FastaReader("/Users/franzweisel/Documents/Systembio/Cutibacterium_granulosum_NCTC11865/GCA_900186975.1_50569_F01_cds_from_genomic.fna")