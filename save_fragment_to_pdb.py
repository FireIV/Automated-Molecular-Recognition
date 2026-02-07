import os
import sys
from Bio.PDB import PDBParser, NeighborSearch, PDBIO
from rdkit import Chem

def parse_structure(pdb_path):
    """Parse a PDB file with BioPython."""
    parser = PDBParser(QUIET=True)
    structure_id = os.path.basename(pdb_path).split('.')[0]
    structure = parser.get_structure(structure_id, pdb_path)
    return structure

input_path = r"C:\Users\bltit\OneDrive\Desktop\Fentanyl Project\4XP6_binding_residues\.pdb Files\4XP6 (2).pdb"
structure = parse_structure(input_path)

model = structure[0]
chain = model["A"]
residue = chain[43]

print(residue.get_full_id())
print(residue.get_resname())
print(list(residue.get_atoms()))



### GETTING THE FILE WHERE IT SHOULD GO
filename = "testoutput"
io = PDBIO()
io.set_structure(residue)
io.save(fr"C:\Users\bltit\OneDrive\Desktop\Molecular Modeling\{filename}.pdb")