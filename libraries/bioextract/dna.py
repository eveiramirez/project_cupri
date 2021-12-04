"""
MODULE NAME
        dna

VERSION
        [0.0.1]

PYTHON VERSION
        3.9

AUTHORS
        Diego

CONTACT
        iramirez@lcg.unam.mx

DESCRIPTION
        Modulo para la obtencion de informacion de genomas.

CATEGORY
        DNA
"""

from bs4 import BeautifulSoup
import requests
import json

failed = []
main = {}
with open('archaea.txt') as m:
    main = json.load(m)


def find_info(key):
    if main[key]['link'][0:4] != "http":
        main[key]['link'] = "https://www.ncbi.nlm.nih.gov/Taxonomy" \
                            "/Browser/" + main[key]['link']
    r = requests.get(main[key]['link'])

    main_page = BeautifulSoup(r.content, 'html5lib')

    information = main_page.find("form").find("table").find_next_sibling().find_next_sibling().find_next_sibling().find_next_sibling().find("td")

    lin = information.find("dd").find_all("a")
    lineage = []
    for taxa in lin:
        lineage.append(taxa.text)

    rank = information.find("h2").find_next_sibling().find_next_sibling().find_next_sibling().find_next_sibling().find_next_sibling().find_next_sibling().text

    main[key]["lineage"] = lineage
    main[key]["rank"] = rank

    bacdive = main_page.find(href="http://bacdive.dsmz.de/")
    gold = main_page.find(href="http://genomesonline.org")
    oma = main_page.find(href="http://omabrowser.org")
    print("general success!")
    try:
        bacdive_link = bacdive.parent.parent.find("a")["href"]
        bd = requests.get(bacdive_link)
        bacdive_page = BeautifulSoup(bd.content, 'html5lib')
        sections = bacdive_page.find(id="content")
        # section 1 (no subtables needed)
        oxygen_tolerance = sections.find(attrs={"data-src-tbl": "oxygen_tolerance"})
        if oxygen_tolerance != None:
            oxygen_tolerance = oxygen_tolerance.find_next_sibling().find_next_sibling().text

        # section 2
        table = sections.find(id="temp_table") #subtable1
        culture_temperature = table.find(attrs={"data-src-tbl": "culture_temp"})
        if culture_temperature != None:
            culture_temperature = culture_temperature.find_next_sibling().find_next_sibling().find_next_sibling().find_next_sibling().text
        table = sections.find(id="temp_range") #subtable2
        temperature_range = table.find(attrs={"data-src-tbl": "culture_temp"})
        if temperature_range != None:
            temperature_range = temperature_range.find_next_sibling().find_next_sibling().text

        # section 3
        table = sections.find(attrs={"class": "id_4 section expandsection"}).find_all("tr")
        location = table[1].find(attrs={"data-src-tbl": "origin"})
        if location != None:
            location = location.find_next_sibling().find_next_sibling().text
        country = table[2].find(attrs={"data-src-tbl": "origin"})
        if country != None:
            country = country.find_next_sibling().find_next_sibling().text
        continent = table[3].find(attrs={"data-src-tbl": "origin"})
        if continent != None:
            continent = continent.find_next_sibling().find_next_sibling().text

        # section 4
        biosafety_level_A = sections.find(attrs={"data-src-tbl": "risk_assessment"}).find_next_sibling().find_next_sibling().text
        """
        print(oxygen_tolerance)
        print(culture_temperature)
        print(temperature_range)
        print(location)
        print(country)
        print(continent)
        print(biosafety_level_A)
        """
        main[key]["oxygen"] = oxygen_tolerance
        main[key]["temperature"] = culture_temperature
        main[key]["temperature_range"] = temperature_range
        main[key]["location_1"] = location
        main[key]["location_2"] = country
        main[key]["location_3"] = continent
        main[key]["biosafety_A"] = biosafety_level_A
        print("b success!")
    except:
        print("failed to retrieve bacdive data for: " + k + ":" +
              main[k]["name"])
        failed.append("failed to retrieve bacdive data for: " + k + ":"
                      + main[k]["name"])
    try:
        gold_link = gold.parent.parent.find("a")["href"]
        gd = requests.get(gold_link)
        gold_page = BeautifulSoup(gd.content, 'html5lib')
        table = gold_page.find(id="OrganismInformation").find("tbody").find_all("tr")

        organism_type = table[18].find("td").find_next_sibling().text.strip()
        is_cultured = table[19].find("td").find_next_sibling().text.strip()
        culture_type = table[20].find("td").find_next_sibling().text.strip()
        biosafety_level_b = table[22].find("td").find_next_sibling().text.strip()
        is_public = table[27].find("td").find_next_sibling().text.strip()
        """
        print(organism_type)
        print(is_cultured)
        print(culture_type)
        print(biosafety_level_b)
        print(is_public)
        """
        main[key]["organism_type"] = organism_type
        main[key]["is_cultured"] = is_cultured
        main[key]["culture_type"] = culture_type
        main[key]["biosafety_B"] = biosafety_level_b
        main[key]["is_public"] = is_public
        print("g success!")
    except:
        print("failed to retrieve gold data for: " + k + ":" +
              main[k]["name"])
        failed.append("failed to retrieve gold data for: " + k + ":" +
                      main[k]["name"])
    try:
        oma_link = oma.parent.parent.find("a")["href"]
        oa = requests.get(oma_link)
        oma_page = BeautifulSoup(oa.content, 'html5lib')
        body = oma_page.find_all("dd")

        sequence_number = body[-4].string
        matrix_protein = body[-3].string
        amino = body[-2].string
        protein_length = body[-1].string
        """
        print(sequence_number)
        print(matrix_protein)
        print(amino)
        print(protein_length)
        """
        main[key]["sequence_number"] = sequence_number
        main[key]["matrix_protein"] = matrix_protein
        main[key]["amino"] = amino
        main[key]["protein_length"] = protein_length
        print("o success!")
    except:
        print("failed to retrieve oma data for: " + k + ":" +
              main[k]["name"])
        failed.append("failed to retrieve oma data for: " + k + ":" +
                      main[k]["name"])


count = 1
for k in main:
    try:
        find_info(k)
        print(count)
    except:
        print("failed to retrieve general data for: " + k)
        failed.append("failed to retrieve general data for: " + k)
    count += 1

with open('archaea_master_dict.txt', 'w') as n:
    json.dump(main, n)

with open('failed_data.txt', 'w') as n:
    json.dump(failed, n)
