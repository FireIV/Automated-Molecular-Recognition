from rcsbapi.search import AttributeQuery
def find_ligand_complexes_rcsb(target_ligand):
    ### Takes ligand and searches rcsb.org for files that contain it ###
    q = AttributeQuery(
        attribute="rcsb_id",
        operator="exact_match",
        value=target_ligand,
        service="text_chem"
    )
    return list(q())

import requests
def download_pdbs(pdb_list):
    ### Takes a list of pdb codes and downloads their respective .pdb files ###
    for pdb in pdb_list:
        file = requests.get(f"https://files.rcsb.org/download/{pdb}.pdb")
        with open(f"{pdb}.pdb", "w") as f:
            f.write(file.text)

def main():
    lig_search_inp = input("Input target ligand's 3-letter code (ex: meth -> B40, fentanyl -> 7V7): ")
    print("Working...")
    pdb_list = find_ligand_complexes_rcsb(lig_search_inp)
    print(f"Found {len(pdb_list)} files!")
    for file in pdb_list:
        idx = 1
        print(f"[{idx}]: {file}")
        idx = idx + 1
    pull_file_req = input("Download these files? [Y/N]: ")
    if pull_file_req == "Y":
        download_pdbs(pdb_list)
    



pdb_id = "4HHB"
url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
response = requests.get(url)

with open(f"{pdb_id}.pdb", "w") as f:
    f.write(response.text)

if __name__ == "__main__":
    main()