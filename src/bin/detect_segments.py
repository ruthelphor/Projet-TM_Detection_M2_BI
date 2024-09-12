#!/usr/bin/python3

"""Module to calculate the mebrane plane on the protein."""

# import modules
import pandas as pd
import math
import matplotlib.pyplot as plt
import os
import shutil


def fibonacci_sphere(samples=10):
    """Generate uniforme distribution of points on a hemisphere.

    This function was taken from the internet.

    Parameters
    ----------
    samples : int
        The number of points wanted on the hemisphere

    Returns
    -------
    list of list
        The coordinates of each point on the hemisphere
    """
    points = []
    phi = math.pi * (3. - math.sqrt(5.))  # golden angle in radians

    for i in range(samples):
        coord_y = 1 - (i / float(samples - 1))  # y goes from 1 to 0
        radius = math.sqrt(1 - coord_y * coord_y)  # radius at y

        theta = phi * i  # golden angle increment

        coord_x = math.cos(theta) * radius
        coord_z = math.sin(theta) * radius

        points.append([coord_x, coord_y, coord_z])

    return points


def position(sphere_points, ca_data):
    """Calculate the hydrophobicity within a slice.

    Parameters
    ----------
    sphere_points : list of list
        The coordinates of points on a hemisphere
    ca_data : Pandas DataFrame
        The coordinates of alpha carbon of a protein and their residu's name

    The DataFrame also contains the accessibility of each atom.

    Returns
    -------
    int
        The best score of hydrophobity
    """
    membrane = {"score": 0,"point": 0, "low_memb":0}
    # Loop on each point of the sphere
    for point in sphere_points:
        mini = pd.Series([ca_data["x"].min(), ca_data["y"].min(), \
                    ca_data["z"].min()], index=["x","y", "z"])  # minimal point
        maxi = pd.Series([ca_data["x"].max(), ca_data["y"].max(), \
                    ca_data["z"].max()], index=["x","y", "z"])  # maximal point
        # loop on each position of the membrane
        while mini["y"] < maxi["y"]:
            nb_ca, hydro = 0, 0
            for i in range(ca_data.shape[0]):  # loop on the CA
                if is_in_membrane(point, ca_data.iloc[i], mini, 15):
                    nb_ca += 1
                    if is_hydrophobe(ca_data.iloc[i]):
                        hydro += 1
            if nb_ca != 0:
                score  = hydro / nb_ca
                if membrane["score"] < score:
                    membrane["score"] = score
                    membrane["point"] = point
                    membrane["low_memb"] = mini
            mini += 1
    return membrane


def is_hydrophobe(carbon):
    """Test if an atom is in a hydrophobe residu.

    Parameters
    ----------
    carbon : Pandas Series
        Series containing the name of its residu

        The Series contain also the coordinate of the atom
        and its solvant accessobility

    Returns
    -------
    Boolean
        If the atom is or not in a hydrophobe residu.
    """
    hydrophobe = ["TRP", "ISO", "LEU", "PHE", "ALA", "MET", "VAL"]

    if carbon["resid_name"] in hydrophobe:
        return True
    return False


def is_in_membrane(sphere_pt, ca_coords, point_min, memb_width):
    """Test if an atom is between two planes.

    Parameters
    ----------
    sphere_pt : list
        The coordinate of a point on the hemisphere
    ca_coords : Pandas Series
        The coordinate of an atom
    point_min : Pandas Series
        The coordinate of the lowest point
    memb_width : int
        The width of the membrane

    Returns
    -------
    Boolean
        If the atom is or not between the planes.
    """
    plane1_d = calculate_plane(sphere_pt, point_min, 0)
    plane2_d = calculate_plane(sphere_pt, point_min, memb_width)
    ca_d = calculate_plane(sphere_pt, ca_coords, 0)

    if plane1_d < ca_d < plane2_d:
        return True
    return False


def calculate_plane(sphere_pt, atom, memb_width):
    """Calculate the plane perpendicular to the sphere point's axis.

    Parameters
    ----------
    sphere_point : list
        The coordinate of a point on the hemisphere
    atom : Pandas Series
        The coordinate of an atom
    memb_width : int
        The width of the membrane

    Returns
    -------
    list of tuple
        The coordinates of each point on the hemisphere
    """
    plane_d = sphere_pt[0] * (atom["x"] + sphere_pt[0] * memb_width) + \
            sphere_pt[1] * (atom["y"] + sphere_pt[1] * memb_width) + \
            sphere_pt[2] * (atom["z"] + sphere_pt[2] * memb_width)
    return plane_d


def plot(ca_data, memb):
    """Plot the protein and the membrane plane.
    
    Stores the results on a png file.

    Parameters
    ----------
    ca_data : Pandas DataFrame
        The coordinates of the alpha carbon in the protein
    memb : dictionnary
        Dictionnary containing the score of the membrane plane,
        the position of the lower and upper layer.
    """
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    
    ax.scatter3D(ca_data["x"], ca_data["y"], ca_data["z"])
    ax.scatter3D(memb["low_memb"]["x"], memb["low_memb"]["y"], 
                memb["low_memb"]["z"])
    ax.scatter3D(memb["low_memb"]["x"] + memb["point"][0]*15, 
                memb["low_memb"]["y"] + memb["point"][1]*15, 
                memb["low_memb"]["z"] + memb["point"][2]*15)

    plt.savefig("protein_with_memb.png", format="png")

    # Obtenir le chemin absolu du répertoire parent
    parent_dir = os.path.abspath("/home/etudiant/Projet_TMDET_reprograming/")

    # Créer le chemin vers le répertoire "resultat"
    resultat_dir = os.path.join(parent_dir, "results")

    # Chemin d'origine du fichier généré dans le répertoire actuel
    fichier_a_deplacer = os.path.join(os.getcwd(), "protein_with_memb.png")

    #  Définir le nouveau chemin où le fichier sera déplacé
    chemin_destination = os.path.join(resultat_dir, "protein_with_memb.png")

    # Déplacer le fichier généré vers le répertoire "resultat"
    shutil.move(fichier_a_deplacer, chemin_destination)

def save_ca_and_membrane_to_pdb(ca_data, membrane, output_pdb):
    """Sauvegarde les atomes alpha carbone et le plan de la membrane sous forme de disques dans un fichier PDB.

    Parameters
    ----------
    ca_data : Pandas DataFrame
        Données des atomes alpha carbone accessibles, incluant les coordonnées.
    membrane : dict
        Contient le score de la membrane, la direction du plan et les coordonnées des points du plan.
    output_pdb : str
        Le nom du fichier PDB de sortie.
    """
    with open(output_pdb, 'w') as f:
        f.write("HEADER    Accessible Alpha Carbons and Membrane Plane\n")

        # Sauvegarder les atomes alpha carbone
        for i, row in ca_data.iterrows():
            f.write(f"ATOM  {i+1:5d}  CA  ALA A{i+1:4d}    {row['x']:8.3f}{row['y']:8.3f}{row['z']:8.3f}  1.00  0.00\n")

        # Ajouter des pseudoatomes pour représenter les bords du plan de la membrane sous forme de disques
        low_point = membrane["low_memb"]
        direction = membrane["point"]

        # Générer un cercle de points autour de low_point et high_point
        num_points = 30  # Nombre de points pour former le disque
        radius = 10  # Rayon du disque

        # Ajouter des points pour le plan bas de la membrane (cercle autour de low_point)
        for j in range(num_points):
            angle = 2 * 3.14159 * j / num_points
            x_disk = low_point["x"] + radius * direction[0] * (j / num_points) * 15
            y_disk = low_point["y"] + radius * direction[1] * (j / num_points) * 15
            z_disk = low_point["z"] + radius * direction[2] * (j / num_points) * 15

            f.write(f"HETATM{len(ca_data) + j + 1:5d}  DSK  MEM B{len(ca_data) + j + 1:4d}    {x_disk:8.3f}{y_disk:8.3f}{z_disk:8.3f}  1.00  0.00\n")

        # Ajouter des points pour le plan haut de la membrane (cercle autour de high_point)
        high_point = {
            "x": low_point["x"] + direction[0] * 15,
            "y": low_point["y"] + direction[1] * 15,
            "z": low_point["z"] + direction[2] * 15
        }

        for j in range(num_points):
            angle = 2 * 3.14159 * j / num_points
            x_disk = high_point["x"] + radius * direction[0] * (j / num_points) * 15
            y_disk = high_point["y"] + radius * direction[1] * (j / num_points) * 15
            z_disk = high_point["z"] + radius * direction[2] * (j / num_points) * 15

            f.write(f"HETATM{len(ca_data) + num_points + j + 1:5d}  DSK  MEM C{len(ca_data) + num_points + j + 1:4d}    {x_disk:8.3f}{y_disk:8.3f}{z_disk:8.3f}  1.00  0.00\n")

        # Terminer le fichier PDB
        f.write("END\n")
