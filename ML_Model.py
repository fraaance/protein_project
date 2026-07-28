# ML - Model
with open("/Users/franzweisel/Documents/Systembio/Cutibacterium_granulosum_TM11/cuti_TM11_translated_cds.faa") as f:
    counter = sum(1 for l in f if l.startswith(">"))
    print(counter)

