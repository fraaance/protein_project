from dataclasses import dataclass
import re


class FastaReader: 
    def __init__(self, file_path):
        self.path = file_path
        input = self.validate_path(file_path)
        self.print_out(input)

    def validate_path(self, file_path):
        if (file_path.endswith(".faa")):
            print("*" * 50)
            print(" FASTA - Amino Acid ".center(50, '*'))
            print("*" * 50)

            header, seq = self.read_file(file_path)
            accession, gene_name = re.split(r"\s+", header[0])
            protein = Protein(accession=accession[1:], 
                              gene_name=gene_name, 
                              organism=header[1].replace("organism=", ""), 
                              gene_id=int(header[2].replace("GeneID=", "").rstrip("]")), 
                              sequence=seq)
            
            return protein
        
        elif (file_path.endswith(".fna")):
            print("*" * 50)
            print(" FASTA - mRNA ".center(50, '*'))
            print("*" * 50)

            header, seq = self.read_file(file_path)
            
            #print(re.split(r"[:\-\s]+|(?<=c)(?=\d)", header[0]))
            accession, strand, start, end, gene_name = re.split(r"[:\-\s]+|(?<=c)(?=\d)", header[0])
            gene_id=int(header[2].replace("GeneID=", "").rstrip("]"))
            rna = RNA(accession=accession[1:], 
                      start=start, end=end, 
                      strand=strand, 
                      gene_name=gene_name, 
                      organism=header[1].replace("organism=", "").rstrip("]"), 
                      gene_id=gene_id, 
                      sequence=seq)
            print(f"\nsequence {gene_name} [GeneID={gene_id}] successfully read".center(50, ' '))
            return rna
        
        else:
            print(" No Valid File - Exit ".center(50, '*'))
            return SystemExit("no valid file!"), False
    
    def read_file(self, file_path):
        with open(file_path) as f:
            lines = f.read().splitlines()
            #header = lines[0]
            header = re.split(r"\s\[+", lines[0])
            #print(len(header))
            seq = "".join(lines[1:])
        return header, seq

    def print_out(self, input):
        print("*" * 50)
        print(input.accession, input.gene_name, input.gene_id, input.organism)

        
@dataclass
class RNA:
    accession: str
    start: int
    end: int 
    strand: str
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