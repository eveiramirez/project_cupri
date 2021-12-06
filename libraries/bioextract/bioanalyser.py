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
from argparse import ArgumentParser, RawTextHelpFormatter
from re import match, search
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter


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
                          download_dir=None):
    """Obtener los reportes de estadisticas de los organismos de la
    base de datos de Assembly. En caso de crearse archivos, se les
    asignara el nombre de acuerdo al ID obtenido de Assembly"""

    # Definir el correo necesario para la busqueda
    Entrez.email = email

    # Inicializar el handle
    handle = None

    # Realizar la busqueda de ids
    ids_orgs = get_ids(Entrez.email, terms, db="assembly")

    # Obtener los URLs de los reportes de estadisticas
    stats_url_list = []

    for id_org in ids_orgs:
        # Obtener el resumen del record
        handle = Entrez.efetch(db="assembly", id=id_org,
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


def get_stats_dictionary(reports: list):
    """Obtener a partir de una lista de reportes de estadisticas de
    ensambles obtenidos por assembly_stats_report una lista de
    diccionarios que almacenen las estadisticas correspondientes a
    todas las secuencias y del reporte"""

    # Creas lista de diccionarios
    stats_list = []

    # Obtener elementos de Assembly Statistics Report
    for organism in reports:
        # Crear diccionario
        stats_dict = {}

        # Obtener los estadisticos antes de Assembly-Units
        last_index = organism.index("## Assembly-Units:")

        for stat in organism[1:last_index - 1]:
            # Separar el nombre del estadistico y su valor
            elements = stat.split(":")

            # Obtener el nombre del estadistico
            elements[0] = elements[0].replace("# ", "")

            # Eliminar espacios que esten en extremos del valor
            elements[1] = elements[1].lstrip()

            # Guardar en el diccionario
            stats_dict[f"{elements[0]}"] = elements[1]

        # Obtener el inicio de las estadisticas de todas las secuencias
        first_index = organism.index("# unit-name\tmolecule-name\t"
                                     "molecule-type/loc\tsequence-type"
                                     "\tstatistic\tvalue")

        # Obtener la posicion de la estadistica final
        last_index = last_match(r'^(all\t){4}', organism)

        # Guardar las estadisticas de todas las secuencias
        for stat in organism[first_index + 1:last_index + 1]:
            stat = stat.replace('all\tall\tall\tall\t', "")

            # Separar el nombre del estadistico y su valor
            elements = stat.split("\t")

            # Guardar en el diccionario el estadistico
            stats_dict[f"{elements[0]}"] = elements[1]

        # Guardar el diccionario en la lista
        stats_list.append(stats_dict)

    return stats_list


def last_match(pattern, list_strings: list[str]):
    """Obtener de una lista la posicion del ultimo elemento que
    coincida con un patron"""

    for stat in list_strings[::-1]:
        if match(pattern, stat):
            return list_strings.index(stat)


def stats_dataframe(email: str, terms: list[str], output=None):
    """Crea un dataframe de las estadisticas de ensambles de
    organismos a partir de una lista con diccionarios"""

    # Obtener las estadisticas
    stats_reports = assembly_stats_report(email, terms)

    # Obtener los diccionarios de las estadisticas
    dictionaries = get_stats_dictionary(stats_reports)

    # Crear la lista de las filas
    i_list = []

    # Obtener todos los indices de las filas
    for stats in dictionaries:
        i_list = list(set(list(stats.keys())+i_list))

    # Guardar las estadisticas de los ensambles
    stats_df = pd.DataFrame(index=i_list)

    for stats in dictionaries:
        stats_series = pd.Series(stats.values(), index=stats.keys(),
                                 name=stats["Assembly name"])

        stats_df[stats_series.name] = stats_series

    # Evaluar si se se guardara el dataframe en un archivo csv

    if type(output) == str:
        try:
            # Crear el archivo del dataframe en formato csv
            stats_df.to_csv(output)
        # Imprimir un error si la ruta del output es invalida
        except IOError as ex:
            print('El archivo no pudo ser creado: ' + ex.strerror)

    else:
        return stats_df


def stats_graph(email: str, terms: list[str], stat: str, output):
    """Crea una grafica de una de las estadisticas de los ensambles de
        organismos"""

    # Obtener el dtaframe con los valores de cada ensamble
    stats_df = stats_dataframe(email, terms)

    # Obtener los valores del estadistico
    values = stats_df.loc[stat].tolist()

    # Verificar si se trata de un valor numerico
    is_num = 0
    for value in values:
        if is_num == 0:
            if type(value) != float:
                if value.isdigit() and stat != "Taxid":
                    is_num = 1

    # Evaluar si se trata de un valor numerico
    if is_num:
        # Obtener los nombres de ensambles
        assemblies_names = stats_df.loc[stat].index.values

        # Cambiar los valores de tipo string a int
        values = list(map(int, values))

        # Generar la figura
        fig, ax = plt.subplots()

        # Generar las barras
        ax.bar(np.arange(len(assemblies_names)), values)

        # Generar los nombres de las barras
        ax.set_xticks(np.arange(len(assemblies_names)),
                      labels=assemblies_names)

        # Añadir etiquetas de los ejes y titulo
        ax.set(xlabel="assemblies", title=stats_df.loc[
            stat].name)

    # Generar grafica del numero de veces que aparece cada valor
    else:
        # Cambiar los valores de tipo string a int
        val_occurrences = Counter(values)

        stats_vals = list(val_occurrences.keys())
        val_occurrences = list(val_occurrences.values())

        # Generar la figura
        fig, ax = plt.subplots()

        # Generar las barras
        ax.bar(np.arange(len(stats_vals)),
               val_occurrences)

        # Generar los nombres de las barras
        ax.set_xticks(np.arange(len(stats_vals)),
                      labels=stats_vals)

        # Añadir etiquetas de los ejes y titulo
        ax.set(xlabel=stats_df.loc[
            stat].name, ylabel="no. of occurrences",
               title=f"no. of occurrences: {stats_df.loc[stat].name}")

    # Rotar las etiquetas de las barras
    plt.setp(ax.get_xticklabels(), rotation=30,
             horizontalalignment='right')
    # Ajustar la imagen
    fig.tight_layout()

    try:
        # Guardar figura
        fig.savefig(output)

    # Imprimir un error si la ruta del output es invalida
    except IOError as ex:
        print('El archivo no pudo ser creado: ' + ex.strerror)


class AssemblyStatisticsReport(object):
    """
    Clase para almacenar el reporte de las estadisticas de un ensamble.

    Atributos
    ---------
    stats: Diccionario de estadisticas del ensamble

    Metodos
    -------
    table(self)
        Genera una tabla con las estadisticas

    """

    def __init__(self, stats: dict):
        self.stats = stats
        self.assembly_name = stats["Assembly name"]
        self.organism_name = stats["Organism name"]

    def table(self):
        """Genera una tabla con las estadisticas"""
        stats_series = pd.Series(self.stats.values(),
                                 index=self.stats.keys(),
                                 name=self.assembly_name)
        return stats_series


# Crear la excepcion AmbiguousBaseError
class AmbiguousError(Exception):
    """Clase para el manejo de errores de distintos tipos."""
    pass


# Crear la descripcion del programa y asignar clase de formato para
# poder controlar la descripcion de ayuda
parser = ArgumentParser(description="Programa que cuenta con "
                                    "diferentes funciones "
                                    "para el analisis y "
                                    "extraccion de datos de "
                                    "NCBI",
                        formatter_class=RawTextHelpFormatter)

# Anadir los argumentos
parser.add_argument("-f", "--function",
                    type=str,
                    help="Funcion a realizar. Requiere que se indique "
                         "el numero de la funcion.\n"
                         "0: get_ids\n"
                         "1: get_tax_data\n"
                         "2: assembly_stats_report\n"
                         "3: stats_dataframe\n"
                         "4: stats_graph")

parser.add_argument("-e", "--email",
                    type=str,
                    help="email requerido para el uso de los servers "
                         "de E-utility")

parser.add_argument("-g", "--organisms",
                    type=str,
                    help="Lista de terminos con los cuales se "
                         "realizara la "
                         "busqueda separados por comas.\nSolo se "
                         "requiere para las "
                         "funciones 0,2,3,4")

parser.add_argument("-i", "--ids",
                    type=str,
                    help="Lista de Taxonomy IDs separados por comas de "
                         "los organismos a buscar.\nEste parametro "
                         "solo es necesario para la funcion 1")

parser.add_argument("-d", "--data_base",
                    type=str,
                    help="Base de datos a realizar la busqueda. Solo "
                         "es necesario para la funcion 0")

parser.add_argument("-o", "--output",
                    type=str,
                    help="Directorio para guardar los archivos "
                         "generados por las funciones 2,3,4. \nEn caso "
                         "de no proporcionarse la funcion 2 devolvera"
                         " los resultados en forma de\nlista, y la "
                         "funcion"
                         " 3 en un dataframe. Para la funcion 4 es "
                         "obligatorio")

parser.add_argument("-s", "--stat",
                    type=str,
                    help="Estadistico a buscar, requerido para la "
                         "funcion 4")

# Ejecutar el metodo parse_args()
args = parser.parse_args()

# Comprobar si se tienen todos los argumentos necesarios
val_errs = []
try:
    argument = [args.function, args.email, args.organisms, args.ids,
                args.data_base, args.output, args.stat]
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
    if argument[0] == "0" or argument[0] == "2" or argument[0] == "3"\
            or argument[0] == "4":
        if argument[2] is None:
            val_errs.append("No se obtuvo la lista de organismos")
            invalid = 1
    if argument[0] == "4":
        if argument[5] is None:
            val_errs.append("No se obtuvo la direccion para la imagen")
            invalid = 1
        if argument[6] is None:
            val_errs.append("No se obtuvo el estadistico a buscar")
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
    # Evaluar si se generara un archivo
    if args.output is None:
        print(assembly_stats_report(args.email, args.organisms.split(
            ",")))
    else:
        assembly_stats_report(args.email, args.organisms.split(
            ","), args.output)
elif function == 3:
    # Evaluar si se generara un archivo
    if args.output is None:
        print(stats_dataframe(args.email, args.organisms.split(
            ",")))
    else:
        # Revisar que el archivo tenga la extension csv
        if search(r".csv$", args.output):
            stats_dataframe(args.email, args.organisms.split(
                ","), args.output)
        else:
            raise AmbiguousError("La extension del archivo no es "
                                 "csv")
elif function == 4:
    stats_graph(args.email, args.organisms.split(
                ","), args.stat, args.output)
else:
    raise SystemExit(f"El numero de funcion '{function}' no pertenece a"
                     f" ninguna funcion")
