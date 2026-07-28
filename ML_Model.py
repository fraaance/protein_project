# ML - Model
import re
with open("/Users/franzweisel/Documents/Systembio/Cutibacterium_granulosum_TM11/cuti_TM11_translated_cds.faa") as f:
    liste = re.split(">", f.read())[1:]
    print(liste[0])

    for e in liste:
        lines = e.splitlines()
        header = lines[0]
        seq = "".join(lines[1:])
        print("header: ", header)
        #print(re.split(header))
    

