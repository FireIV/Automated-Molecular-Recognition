from rdkit import Chem
import sys
pdb_path = input(f"Paste filepath for .pdb file here: ").strip("\"")
rdkit_structure = Chem.rdmolfiles.MolFromPDBFile(pdb_path,True,True,0,True)
conf = rdkit_structure.GetConformer()
for atom in rdkit_structure.GetAtoms():
    pos = conf.GetAtomPosition(atom.GetIdx())
    print(f"{atom.GetSymbol()}{atom.GetIdx()}: {pos.x}, {pos.y}, {pos.z}")

add_hs = ""
while add_hs != "y":
    add_hs = input("Would you like to add hydrogens? [y/n]: ").lower()
    if add_hs == "n":
        print("Program Complete!")
        sys.exit()
    if add_hs != "y":
        add_hs = input("Invalid option, please choose if you'd like to add hydrogens. [y/n]: ")


# If you remove this block (lines 21-40) the code functions as designed, this is just unfinished
acid_base_zwit_check = ""
change_valency = ""
while acid_base_zwit_check != "n":
    acid_base_zwit_check = input("Are there any nonstandard valences in this molecule? [y/n]: ").lower
    if acid_base_zwit_check == "n":
        break
    if acid_base_zwit_check != "y":
        print("Invalid Input")
    if acid_base_zwit_check == "y":
        for atom in rdkit_structure.GetAtoms():
            print(f"[{atom.GetIdx}] {atom.GetSymbol()}")
        while change_valency != "n":
            change_atom_valency = input(f"Please select the index number of the atom whose formal charge you'd like to change [0-{len(atom in rdkit_structure.GetAtoms())}], or if you'd like to exit this mode, type [n]: ")
            if change_atom_valency == "n":
                break
            try:
                int(change_atom_valency)
            except:
                print("Invalid input")
            formal_charge = input("What is the desired formal charge for this atom? :")
            ### Here I want to be able to change the formal charge of atoms for nonstandard valencies, like methamphetamine being a conjugate acid
            ### C(NH_2^+)C in vivo, and I planned on doing this through atom.SetFormalCharge(int) but as you can see... this seems like a bit of a 
            ### mess with all the while loops and such, and so it's as far as I've gotten.

        
with_hydrogens = Chem.rdmolops.AddHs(rdkit_structure, addCoords=True, addResidueInfo=True)
conf_withH = with_hydrogens.GetConformer()
for atom in with_hydrogens.GetAtoms():
    pos = conf_withH.GetAtomPosition(atom.GetIdx())
    print(f"{atom.GetSymbol()}{atom.GetIdx()}: {pos.x}, {pos.y}, {pos.z}")

Chem.MolToPDBFile(with_hydrogens,'testfile', flavor = 2)