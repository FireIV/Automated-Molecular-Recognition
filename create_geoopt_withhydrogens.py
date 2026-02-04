import os
import sys
from Bio.PDB import PDBParser, NeighborSearch, PDBIO
from rdkit import Chem

#Get input file from user
def parse_structure(pdb_path):
    """Parse a PDB file with BioPython."""
    parser = PDBParser(QUIET=True)
    structure_id = os.path.basename(pdb_path).split('.')[0]
    structure = parser.get_structure(structure_id, pdb_path)
    return structure

#Identify target ligand
def find_hetero_residues(structure):
    """
    Find all hetero residues in the structure except water.
    BioPython residue.id format: (hetflag, resseq, icode)
    """
    hetero_residues = []
    for model in structure:
        for chain in model:
            for residue in chain:
                hetflag, resseq, icode = residue.id
                if isinstance(hetflag, str) and len(hetflag) > 0 and hetflag[0] == "H":
                    # Exclude water molecules
                    if residue.get_resname().upper() not in ("HOH", "WAT"):
                        hetero_residues.append((model, chain, residue))
    return hetero_residues

def describe_hetero_residues(hetero_residues):
    """Produce human-readable labels for hetero residues."""
    labels = []
    for i, (model, chain, residue) in enumerate(hetero_residues, start=1):
        hetflag, resseq, icode = residue.id
        resname = residue.get_resname()
        label = f"{i}: {resname} {hetflag}{resseq}{icode} in chain {chain.id}, model {model.id}"
        labels.append(label)
    return labels

def choose_hetero_residue(hetero_residues):
    """Allow user to select a hetero residue from the list."""
    labels = describe_hetero_residues(hetero_residues)
    if not labels:
        print("No hetero residues found.")
        sys.exit(0)

    print("\nHetero residues found:")
    for label in labels:
        print(label)

    while True:
        het_res_choose = input("\nEnter index of hetero residue to analyze (or 'q' to quit): ")
        if het_res_choose.lower() == "q":
            sys.exit(0)
        try:
            idx = int(het_res_choose)
            if 1 <= idx <= len(hetero_residues):
                return hetero_residues[idx - 1]
        except ValueError:
            pass
        print("Invalid selection. Try again.")


#Locate neighbor residues
def find_neighbor_residues(structure, center_residue, cutoff=6.0):
    """
    Find all residues within cutoff Angstroms of the center_residue
    with BioPython's NeighborSearch.
    """
    atoms = list(structure.get_atoms())
    ns = NeighborSearch(atoms)
    center_atoms = list(center_residue.get_atoms())
    neighbor_residues = set()

    for atom in center_atoms:
        neighbors = ns.search(atom.get_coord(), cutoff, level="R")
        for res in neighbors:
            if res is not center_residue:
                neighbor_residues.add(res)

    sorted_neighbors = sorted(
        neighbor_residues,
        key=lambda r: (r.get_parent().id, r.id[1], r.id[2])
    )
    return sorted_neighbors

def describe_residue(residue):
    """Produce readable description of a residue."""
    hetflag, resseq, icode = residue.id
    resname = residue.get_resname()
    chain_id = residue.get_parent().id
    return f"{resname} {hetflag}{resseq}{icode} chain {chain_id}"

#Pull every residue individually and make them a .pdb file
#Add hydrogens to the ligand (request valency changes if necessary)
#Add hydrogens to each residue (consider charged residues!!)
#Make a Gaussian .com file for each ligand-residue pairing


def main():
    input_pdb = input("Please paste the .pdb file path you'd like to use here: ").strip("\"")
    
    structure = parse_structure(input_pdb)
    
    hetero_residues = find_hetero_residues(structure)
    if not hetero_residues:
        print("No hetero residues detected.")
        return
    
    model, chain, chosen_hetero = choose_hetero_residue(hetero_residues)
    print(f"\nChosen hetero residue: {describe_residue(chosen_hetero)}")
    print("\nSearching for residues within 6 Angstroms...")
    neighbor_residues = find_neighbor_residues(structure, chosen_hetero, cutoff=6.0)
    if not neighbor_residues:
        print("No residues within 6 Å of the chosen hetero residue.")
        return
    print(f"Found {len(neighbor_residues)} neighbor residue(s):")
    for i, res in enumerate(neighbor_residues, start=1):
        print(f"  {i}: {describe_residue(res)}")
    
    plot_frag_chk = input("Would you like to create .pdb files for each ligand/neighbor interaction? [y/n]: ")
    while plot_frag_chk != "y":
        if plot_frag_chk == "n":
            sys.exit()
            plot_frag_chk = input(f"Invalid Input \nWould you like to create .pdb files for each ligand/neighbor interaction? [y/n]: ")
    
    exclude_hetatoms = input("Exclude ligand/heteroatom interactions? [y/n]: ").lower()
    while exclude_hetatoms != "y" and exclude_hetatoms != "n":
        exclude_hetatoms = input(f"Invalid Input \nWould you like to exclude ligand/heteroatom interactions? [y/n]: ")

if __name__ == "__main__":
    main()