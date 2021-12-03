"""
MODULE NAME
        bioanalyser

VERSION
        [0.0.1]

PYTHON VERSION
        3.9

AUTHORS
        Ignacio Emmanuel Ramirez Bernabe
        Diego
        Melissa

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
from urllib.request import (urlopen, urlretrieve)


# Crear la excepcion AmbiguousBaseError
class AmbiguousError(Exception):
    """Clase para el manejo de errores de distintos tipos."""
    pass


def get_ids(email: str, species: list[str], db: str):
    """Obtener los ids de una base de datos"""

    # Definir el correo necesario para la busqueda
    Entrez.email = email

    # Realizar la busquedas de ids
    ids = []
    for specie in species:
        # Guardar la busqueda
        handle = Entrez.esearch(term=specie, db=db,
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


def get_tax_data(email: str, ids: list[str]):
    """Obtener los linajes de los ids de la base de datos de
    Taxonomy. Requiere que la lista de ids sean TaxIDs"""

    # Definir el correo necesario para la busqueda
    Entrez.email = email

    # Obtener los datos de Taxonomy
    handle = Entrez.efetch(db="Taxonomy", id=ids, retmode="xml")

    # Leer el archivo
    organisms = Entrez.read(handle)

    # Obtener los linajes
    lineages = []
    for organism in organisms:
        organism = organism["Lineage"].split(";")
        lineages.append(organism)

    # Devolver los linajes
    return lineages


def assembly_stats_report(email: str, terms: list[str],
                          download_dir=None):
    """Obtener los reportes de estadisticas de los organismos de la
    base de datos de Assembly. En caso de crearse archivos, se les
    asignara el nombre de acuerdo al ID obtenido de Assembly"""

    # Definir el correo necesario para la busqueda
    Entrez.email = email

    # Inicializar handle
    handle = None

    # Realizar la busqueda de ids
    ids_orgs = get_ids(Entrez.email, terms, db="Assembly")
    print(ids_orgs)

    # Obtener los URLs de los reportes de estadisticas
    stats_url_list = []

    for id in ids_orgs:
        handle = Entrez.efetch(db="Assembly", id=id,
                               rettype="docsum")

        organism = Entrez.read(handle)

        # Obtener el url de las estadisticas
        assembly_stats_url = organism['DocumentSummarySet'][
            "DocumentSummary"][0]["FtpPath_Stats_rpt"]
        assembly_stats_url = assembly_stats_url.replace("ftp:",
                                                        "https:")
        stats_url_list.append(assembly_stats_url)

    # Evaluar si se descargaran los archivos o se guardaran en una
    # lista
    if type(download_dir) == str:
        for i, link in enumerate(stats_url_list):
            urlretrieve(link, f"{download_dir}/{ids_orgs[i]}.txt")
    else:
        # Guardar las lineas de cada archivo en un elemento de una lista
        stats_list = []
        for link in stats_url_list:
            # Abrir el archivo y leerlo
            file = urlopen(link)
            file = file.readlines()

            # Cambiar la codificación de las lineas por utf-8,
            # y eliminar los caracteres \r y \n
            for i, line in enumerate(file):
                file[i] = line.decode(encoding='utf-8').replace("\r\n",
                                                                "")
            # Guardar las lineas del archivo en la lista
            stats_list.append(file)

        handle.close()

        return stats_list
