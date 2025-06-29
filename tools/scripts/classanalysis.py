#!/usr/bin/env python3
import re
import sys
import os

#The aim of this script is to browse swig file to get all the included hpp files and then to find
#all the method whose name start with "create" in order to add %newobject instruction in the swig file
#and avoid memory leak (garbage collector will then call delete of the object).

def process_cpp_file(classname,header_file, output_file, first):
    # Catch all functions named 'createFrom*' or 'createSomething*'
    pattern = re.compile(r'\w+\s*\*\s*(create\w+)')
    if first:
        mode = 'w'
    else:
        mode = 'a'
    with open(header_file, 'r', encoding='utf-8') as cpp_file, open(output_file, mode, encoding='utf-8') as output:
        for line in cpp_file:
            match = pattern.search(line)
            if match:
                output.write(f"%newobject {classname}::{match.group(1)};\n")
    return False

def extract_included_files(file_path):
    """
    Parcourt un fichier et extrait tous les noms de fichiers inclus avec #include.

    :param file_path: Chemin vers le fichier à analyser.
    :return: Une liste des noms de fichiers inclus.
    """
    included_files = []
    # Expression régulière pour capturer le contenu entre guillemets après #include
    include_pattern = re.compile(r'%include\s+(\S+\.hpp)\s*$')

    with open(file_path, 'r', encoding='utf-8') as file:
        for line in file:
            match = include_pattern.search(line)
            if match:
                # Ajouter le nom du fichier inclus à la liste
                included_files.append(match.group(1))

    return included_files


if __name__ == "__main__":
    filename  = sys.argv[1]
    output_txt_file = sys.argv[2]
    first = True
    files = extract_included_files(filename)
    for file in files:
        fsplit = file.split("/")
        if len(fsplit) != 2:
            continue
        path = os.path.join(fsplit[0], fsplit[1])
        if path[0] in ["include", "Core"]:
            continue
        classname = fsplit[1].split(".")[0]
        header_file = os.path.join("..","..","include", file)  
        first = process_cpp_file(classname,header_file, output_txt_file, first)
