"""
PACKAGE NAME
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
        bioanalyser es una paqueteria que contiene distintas
        herramientas para el manejo de información de bases de datos
        de NCBI, enfocadas principalmente hacia bacterias

CATEGORY
        Bacteria
        DNA
        Protein
        Sequence

"""

from bioanalyser import (get_ids, get_tax_data, assembly_stats_report,
                         get_stats_dictionary, stats_dataframe,
                         AssemblyStatisticsReport)

__all__ = ['get_ids', 'get_tax_data', 'assembly_stats_report',
           'get_stats_dictionary', 'stats_dataframe',
           'AssemblyStatisticsReport']
