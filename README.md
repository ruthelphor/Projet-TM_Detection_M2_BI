### Projet de Détection des Segments Transmembranaires (TMDET)

# Auteur : Bonou Yves YAMADJAKO

# Description

Ce projet contient un programme conçu pour analyser une protéine fournie sous forme de fichier PDB et identifier les segments transmembranaires. Le programme est une reprogrammation de l'algorithme TMDET pour assigner et détecter les segments transmembranaires, en se basant sur les résidus hydrophobes et les informations structurelles de la protéine (Surface accessible au solvant, hydrophobicité et atomes de carbones alpha). Il génère des fichiers PDB et PNG montrant la structure de la protéine et sa relation avec les plans membranaires.

# Configurer votre environnement

CLoner le dépôt:

git clone https://github.com/CNOV0/Projet_Management.git

Déplacez-vous dans le nouveau répertoire :

cd Projet_TMDET_reprograming

# Prérequis

    Conda (ou Miniconda/Anaconda) pour la gestion des environnements.
    Un fichier PDB de la protéine que vous souhaitez analyser.

# Installation de l'Environnement Conda

Pour installer et configurer l'environnement Conda nécessaire pour exécuter ce projet, utiliser le fichier environment.yml fourni avec ce projet.

# Étapes d'installation de l'environnement :

    Téléchargez le projet et placez-le dans un répertoire de travail.

    Ouvrez un terminal et naviguez vers le répertoire du projet.

    Exécutez la commande suivante pour créer l'environnement Conda à partir du fichier environment.yml :

    conda env create -f environment.yml

# Activation d l'environnement conda
    conda activate Projet_TMP

# Utilisation

    Une fois dans le répertoire du projet, déplacer-vous vers le sous-répertoire src

        cd src/bin/

    Pour exécuter le programme, utiliser la commande suivante:

       python3 main.py <chemin_vers_fichier_PDB> -v

        <chemin_vers_fichier_PDB> : Le chemin vers le fichier PDB de la protéine que vous souhaitez analyser.


        -v : (Optionnel) Activez ce drapeau pour exécuter le programme en mode "verbose" et afficher les étapes détaillées de l'exécution.

    Pour avoir de l'aide, éventuellement à l'utilisation du programme, utiliser la commande

        python3 main.py -h

         ou

        python3 main.py <chemin_vers_fichier_PDB> -h

# Exemple d'exécution

    Voici un exemple d'exécution typique avec tous les paramètres fournis :

        python3 main.py /home/etudiant/Projet_TMDET_reprograming/data/2k8j.pdb -o results/protein_with_memb.png -v

    Dans cet exemple :

        - Le fichier PDB 2k8j.pdb est analysé.
        - Le fichier d'image PNG sera sauvegardé dans le répertoire results/ avec le nom protein_with_memb.png.
        - Le mode "verbose" affichera les détails des étapes dans le terminal.

# Structure des fichiers

    Voici un aperçu des fichiers clés du projet :

    main.py : Script principal qui orchestre l'exécution complète de l'analyse (lecture du fichier PDB, extraction des Cα, calcul des plans membranaires, génération des fichiers de sortie).

    extract_ca_accessible.py : Extraction des carbones alpha accessibles au solvant.

    center_protein.py : Calcul du centre de masse et centrage des coordonnées de la protéine.

    detect_segments.py : Calcul des segments transmembranaires et génération des fichiers de sortie (PNG et PDB).

    environment.yml : Fichier de configuration pour créer l'environnement Conda requis.


        




