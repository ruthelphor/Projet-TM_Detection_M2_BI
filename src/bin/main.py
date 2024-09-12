#!/usr/bin/python3

#Import modules
import sys
import argparse
import pandas as pd
from extract_ca_accessible import find_ca_access, dssp
from center_protein import calculate_center_of_mass, center_protein
from detect_segments import fibonacci_sphere, position, plot, save_ca_and_membrane_to_pdb
import os
import shutil

def main(pdb_file, output_file=None, verbose=False):
    # Function for verbose mode with immediate flush to show steps
    def log(message):
        if verbose:
            print(message, flush=True)  # Force flush for messages to appear immediately
    
    # Run DSSP to get solvent accessibility
    log("[INFO] Running DSSP to get solvent accessibility...")
    dssp_data = dssp(pdb_file)  # Handle DSSP errors within the function
    
    if not dssp_data:
        log("[ERROR] DSSP failed to process the file. Exiting.")
        sys.exit(1)

    # Extract alpha carbons that are accessible to the solvent
    log("[INFO] Extracting solvent-accessible alpha carbons from PDB file...")
    ca_data = find_ca_access(dssp_data, pdb_file)
    log(f"  Number of accessible alpha carbons extracted: {len(ca_data)}")

    # Calculate the center of mass of the protein
    log("[INFO] Calculating the center of mass and centering the protein...")
    center_of_mass = calculate_center_of_mass(ca_data)
    ca_data = center_protein(ca_data, center_of_mass)
    log(f"  Center of mass: {center_of_mass}")

    # Generate directions with the Fibonacci sphere for transmembrane segment detection
    log("[INFO] Generating directions on the sphere...")
    directions = fibonacci_sphere(samples=100)

    # Detect transmembrane segments
    log("[INFO] Detecting transmembrane segments...")
    membrane = position(directions, ca_data)
    log(f"  Best membrane score: {membrane['score']}")

    # Plotting the protein and the membrane plane
    log(f"[INFO] Plotting the protein and membrane to {output_file}...")
    plot(ca_data, membrane)
    log(f"[INFO] Plot saved")

    # Sauvegarder les carbones alpha et la membrane dans un fichier PDB
    output_pdb = "ca_and_membrane.pdb"
    log(f"[INFO] Saving accessible carbons and membrane to {output_pdb}...")
    
    # Appel de la fonction pour sauvegarder les résultats dans un fichier PDB
    save_ca_and_membrane_to_pdb(ca_data, membrane, output_pdb)

    # Obtenir le chemin absolu du répertoire parent
    parent_dir = os.path.abspath("/home/etudiant/Projet_TMDET_reprograming/")

    # Créer le chemin vers le répertoire "resultat"
    resultat_dir = os.path.join(parent_dir, "results")

    # Chemin d'origine du fichier généré dans le répertoire actuel
    fichier_a_deplacer = os.path.join(os.getcwd(), "ca_and_membrane.pdb")

    #  Définir le nouveau chemin où le fichier sera déplacé
    chemin_destination = os.path.join(resultat_dir, "ca_and_membrane.pdb")

    # Déplacer le fichier généré vers le répertoire "resultat"
    shutil.move(fichier_a_deplacer, chemin_destination)


    log("[FINISHED] Program executed successfully.")

if __name__ == "__main__":
    # Manage command-line arguments with help functionality
    parser = argparse.ArgumentParser(
        description="Program for detecting transmembrane segments in a protein.",
        formatter_class=argparse.RawTextHelpFormatter
    )

    # Argument for the PDB file
    parser.add_argument(
        "pdb_file",
        help="The PDB file of the protein to analyze."
    )

    # Optional argument for output file
    parser.add_argument(
        "-o", "--output",
        help="Output file to save the results (plot).",
        default=None
    )

    # Option to activate verbose mode
    parser.add_argument(
        "-v", "--verbose",
        help="Activate verbose mode to show execution details.",
        action="store_true"
    )

    # Parse arguments and run the main program
    args = parser.parse_args()
    main(args.pdb_file, args.output, args.verbose)
