#!/usr/bin/python3

"Module to read a pdb file and extract the alpha carbon that are accessible to the solvant."

# import modules
import pandas as pd
from Bio.PDB import PDBParser
from Bio.PDB import DSSP

def dssp(file_pdb):
    """Run DSSP on a PDB file.

    Used to get the solvant accessibility of each atom.

    Parameters
    ----------
    file_pdb : str
        The name of a PDB file

        The filename may contain de path to fetch the file.

    Returns
    -------
    list ?
        Different informations on the atom
        - especially the accessibility.
    """
    parsed = PDBParser()
    structure = parsed.get_structure("2K1A", file_pdb)
    model = structure[0]
    dssp_data = DSSP(model, file_pdb)

    return dssp_data


def access_solvant(dssp_data, index):
    """Get the value of accessibility for a given atom.

    Parameters
    ----------
    dssp_data : list ?
        A list containing various informations and the accessibility
    index : int
        The atom for which we want the accessibility

    Returns
    -------
    int
        The solvant accessibility of the atom
    """
    a_key = list(dssp_data.keys())[index]
    return dssp_data[a_key][3]  # returns the accessibility score for each CA


def find_ca_access(dssp_data, file_pdb):
    """Find the alpha carbon accessible to the solvant.

    Find the alpha carbon in a PDB file
    and select only those who are accessible to the solvant.

    Parameters
    ----------
    dssp_data : list ?
        A list containing various informations and the accessibility
    file_pdb : str
        The name of a PDB file

        The filename may contain de path to fetch the file.

    Returns
    -------
    Pandas DataFrame
        DataFrame containing the solvant accessibility, the residu's name,
        the coordinates of the alpha carbon.
    """
    carbon_alpha = []

    with open(file_pdb, "r") as file_in:
        for ligne in file_in:
            if ligne.startswith("ATOM") and "CA" in ligne:
                carbon_alpha.append(ligne)
            elif ligne.startswith("ENDMDL"):  # to get just 1 model
                break

    ca_access = []

    for i in range(len(carbon_alpha)):
        access = access_solvant(dssp_data, i)
        if access > 0.5:
            data = [access, carbon_alpha[i].split()[3], \
                    float(carbon_alpha[i].split()[6]), \
                    float(carbon_alpha[i].split()[7]), \
                    float(carbon_alpha[i].split()[8])]  # get z
            ca_access.append(data)

    return pd.DataFrame(ca_access, columns=["accessibility", "resid_name",
                                            "x", "y", "z"])
