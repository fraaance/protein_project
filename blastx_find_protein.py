import subprocess

result = subprocess.run(
    ["blastx", "-query", "query.fasta", "-db", "protein_db"],
    check=True,
    capture_output=True,
    text=True,
)