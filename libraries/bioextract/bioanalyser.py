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
        Modulo principal para analisis de archivos

CATEGORY
        Extractor
        Analyser
        Sequence

"""

from Bio import Entrez
from urllib.request import (urlopen, urlretrieve)
import argparse


# Crear la excepcion AmbiguousBaseError
class AmbiguousError(Exception):
    """Clase para el manejo de errores de distintos tipos."""
    pass


def get_ids(email: str, organisms: list[str], db: str):
    """Obtener los ids de una base de datos"""

    # Definir el correo necesario para la busqueda
    Entrez.email = email

    # Crear la lista para los IDs
    ids = []

    # Realizar la busquedas de ids
    for organism in organisms:
        # Hacer la busqueda y guardarla en el handle
        handle = Entrez.esearch(term=organism, db=db,
                                retmode="xml")
        # Obtener la estructura del archivo XML de la busqueda
        record = Entrez.read(handle)

        # Evaluar si la lista de IDs esta vacia
        if len(record['IdList']) == 0:
            raise SystemExit(
                f"Organismo '{organism}' no encontrado")
        else:
            # Guardar los IDs
            for id_org in record['IdList']:
                ids.append(id_org)

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
                          download_dir):
    """Obtener los reportes de estadisticas de los organismos de la
    base de datos de Assembly. En caso de crearse archivos, se les
    asignara el nombre de acuerdo al ID obtenido de Assembly"""

    # Definir el correo necesario para la busqueda
    Entrez.email = email

    # Inicializar el handle
    handle = None

    # Realizar la busqueda de ids
    ids_orgs = get_ids(Entrez.email, terms, db="Assembly")

    # Obtener los URLs de los reportes de estadisticas
    stats_url_list = []

    for id_org in ids_orgs:
        # Obtener el resumen del record
        handle = Entrez.efetch(db="Assembly", id=id_org,
                               rettype="docsum")

        # Obtener su estructura
        organism = Entrez.read(handle)

        # Obtener el URL de las estadisticas
        assembly_stats_url = organism['DocumentSummarySet'][
            "DocumentSummary"][0]["FtpPath_Stats_rpt"]
        assembly_stats_url = assembly_stats_url.replace("ftp:",
                                                        "https:")
        # Guardar el URL en una lista
        stats_url_list.append(assembly_stats_url)

    handle.close()

    # Evaluar si se descargaran los archivos o se guardaran en una
    # lista
    if type(download_dir) == str:
        # Crear los archivos de los reportes de cada organismo en el
        # directorio indicado
        try:
            for i, link in enumerate(stats_url_list):

                # Guardar el archivo con su UID de Assembly como nombre
                urlretrieve(link, f"{download_dir}/{ids_orgs[i]}.txt")

        # Imprimir un error si la ruta del output es invalida
        except IOError as ex:
            print('Los archivos no pudieron ser creados: ' +
                  ex.strerror)
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

        return stats_list


# Crear la descripcion del programa
parser = argparse.ArgumentParser(description="Programa que cuenta con "
                                             "diferentes funciones "
                                             "para el analisis y "
                                             "extraccion de datos de "
                                             "NCBI")

# Anadir los argumentos
parser.add_argument("-f", "--function",
                    type=str,
                    help="Funcion a realizar. Requiere que se indique "
                         "el numero de la funcion."
                         "0: get_ids"
                         "1: get_tax_data"
                         "2: assembly_stats_report")

parser.add_argument("-e", "--email",
                    type=str,
                    help="email requerido para el uso de los servers "
                         "de E-utility")

parser.add_argument("-g", "--organisms",
                    type=str,
                    help="Lista de terminos con los cuales se "
                         "realizara la "
                         "busqueda separados por comas. Solo se "
                         "requiere para las "
                         "funciones 0,2")

parser.add_argument("-i", "--ids",
                    type=str,
                    help="Lista de Taxonomy IDs separados por comas de "
                         "los organismos a buscar. Este parametro solo "
                         "es necesario para la funcion 1")

parser.add_argument("-d", "--data_base",
                    type=str,
                    help="Base de datos a realizar la busqueda. Solo "
                         "es necesario para la funcion 0")

parser.add_argument("-o", "--output",
                    type=str,
                    help="Directorio para guardar los archivos "
                         "generados por la funcion 2. En caso de no "
                         "proporcionarse la funcion devolvera los "
                         "resultados en forma de lista")

# Ejecutar el metodo parse_args()
args = parser.parse_args()

# Comprobar si se tienen todos los argumentos necesarios
val_errs = []
try:
    argument = [args.function, args.email, args.organisms, args.ids,
                args.data_base, args.output]
    invalid = 0
    if argument[0] is None:
        val_errs.append("No se obtuvo la funcion que se desea ejecutar")
        invalid = 1
    if argument[0] == "0":
        if argument[4] is None:
            val_errs.append("No se obtuvo la base de datos")
            invalid = 1
    if argument[0] == "1":
        if argument[3] is None:
            val_errs.append("No se obtuvo la lista de TaxIDs")
            invalid = 1
    if argument[0] == "0" or argument[0] == "2":
        if argument[2] is None:
            val_errs.append("No se obtuvo la lista de organismos")
            invalid = 1
    if argument[1] is None:
        val_errs.append("No se obtuvo el email")
        invalid = 1
    if invalid:
        raise ValueError
except ValueError:
    for val_err in val_errs:
        print(f"{val_err}: El valor es None")
    raise SystemExit(1)

# Evaluar si para el parametro de funcion se dio un numero valido
arg_str = args.function
if arg_str.isdigit():
    function = int(arg_str)
else:
    raise SystemExit(f"El numero de funcion '{arg_str}' no pertenece a "
                     f"ninguna funcion")


# Evaluar cual funcion ejecutar
if function == 0:
    print(get_ids(args.email, args.organisms.split(","), 
                  args.data_base))
elif function == 1:
    print(get_tax_data(args.email, args.ids.split(",")))
elif function == 2:
    if args.output is None:
        print(assembly_stats_report(args.email, args.organisms.split(
            ","), args.output))
    else:
        assembly_stats_report(args.email, args.organisms.split(
            ","), args.output)
else:
    raise SystemExit(f"El numero de funcion '{function}' no pertenece a"
                     f" ninguna funcion")
