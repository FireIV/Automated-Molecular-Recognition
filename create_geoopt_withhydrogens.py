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

def store_atom_coords(reslist, folder_path):
###Takes a list of reidues and stores each one's atomic coordinates separately
    io = PDBIO()
    for res in reslist:
        io.set_structure(res)
        io.save(rf"{folder_path}{res.get_resname()}{res.get_id()[1]}.pdb")
    
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
    
    plot_frag_chk = input("Would you like to create .pdb files for each ligand/neighbor interaction? [y/n]: ").lower()
    while plot_frag_chk != "y":
        if plot_frag_chk == "n":
            sys.exit()
            plot_frag_chk = input(f"Invalid Input \nWould you like to create .pdb files for each ligand/neighbor interaction? [y/n]: ")
    
    include_hetatms_chk = input("Include other heteroatoms? [y/n]: ").lower()
    while include_hetatms_chk != "y" and include_hetatms_chk != "n":
        include_hetatms_chk = input("Invalid input. Please select yes [y] or no [n]: ")
    
    #need to retool this later to actually remove heteroatoms!!
    #for res in neighbor_residues:
    #    if include_hetatms_chk == "n":
    #       het_chk = (res.get_id()[0])
    #       if het_chk == " ":
    #           print(res.get_full_id()[3])
    #    else:
    #       print(res.get_full_id()[3])
    
    noHs_path = "C:\\Users\\bltit\\Desktop\\test\\"
    store_atom_coords(neighbor_residues, noHs_path)
    def addHs_resfrag(in_path, out_path):
        """
        Add hydrogens to each fragment .pdb file
        
        :param in_path: input folder path of .pdb files without hydrogens
        :param out_path: output folder path of .pdb files with hydrogens added
        """
        for f in os.listdir(in_path):

            f_path = in_path + f
            print(f_path)
            rdkit_struct = Chem.rdmolfiles.MolFromPDBFile(f_path,True,True,0,True)
            print(type(rdkit_struct))
        #    conf = rdkit_struct.GetConformer()
        #    #for atom in rdkit_struct:
        #    #    pos =conf.GetAtomPosition(atom.GetIdx())
        #    #    print(f"{atom.GetSymbol()}{atom.GetIdx()}: {pos.x}, {pos.y}, {pos.z}")
            with_hs = Chem.rdmolops.AddHs(rdkit_struct, addCoords=True, addResidueInfo=True)
            conf_withH = with_hs.GetConformer()
        #    #for atom in with_hs.GetAtoms():
        #    #    pos = conf_withH.GetAtomPosition(atom.GetIdx())
        #    #    print(f"{atom.GetSymbol()}{atom.GetIdx()}: {pos.x}, {pos.y}, {pos.z}")
            f_out = f_path.removesuffix(".pdb") + "_HsAdded.pdb"
            print(f_out)
            Chem.MolToPDBFile(with_hs,f_out, flavor = 2)
    addHs_resfrag(noHs_path,noHs_path)




if __name__ == "__main__":
    main()