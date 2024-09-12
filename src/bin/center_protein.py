#!/usr/bin/python3

#Import modules
import pandas as pd
import sys

def calculate_center_of_mass(ca_data):
    """Calcule le centre de masse des atomes Cα"""
    # Calcul de la moyenne des coordonnées x, y, z
    x_com = ca_data["x"].mean()
    y_com = ca_data["y"].mean()
    z_com = ca_data["z"].mean()
    center_of_mass = [x_com, y_com, z_com]
    return center_of_mass

def center_protein(ca_data, center_of_mass):
    """Centre les coordonnées de la protéine autour de (0, 0, 0)"""
    # Soustraction du centre de masse pour centrer la protéine
    ca_data["x"] -= center_of_mass[0]
    ca_data["y"] -= center_of_mass[1]
    ca_data["z"] -= center_of_mass[2]
    return ca_data 