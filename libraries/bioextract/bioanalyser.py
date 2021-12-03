"""
MODULE NAME
        bioanalyser

VERSION
        [0.0.1]

PYTHON VERSION
        3.9

AUTHORS
        Ignacio Emmanuel Ramirez Bernabe

CONTACT
        iramirez@lcg.unam.mx

DESCRIPTION
        Modulo principal para análisis de archivos

CATEGORY
        Extractor
        Analyser
        Sequence

"""

from Bio import Entrez


# Crear la excepcion AmbiguousBaseError
class AmbiguousError(Exception):
    """Clase para el manejo de errores de distintos tipos."""
    pass


def get_taxid(email: str, species: list[str]):
    """Obtener los TaxID"""

    # Definir el correo necesario para la busqueda
    Entrez.email = email

    # Realizar la busquedas de ids
    ids = []
    for specie in species:
        # Guardar la busqueda
        handle = Entrez.esearch(term=specie, db="taxonomy",
                                retmode="xml")
        # Obtener la estructura del archivo XML de la busqueda
        record = Entrez.read(handle)

        # Evaluar si la lista esta vacia
        if len(record['IdList']) == 0:
            raise SystemExit(
                f"Organismo '{specie}' no encontrado")
        else:
            # Guardar el id
            ids.append(record['IdList'][0])

    # Devolver los IDs
    return ids
