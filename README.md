# project_cupri
https://github.com/eveiramirez/project_cupri

## Introduccion

En este proyecto se trabajo con las herramientas aprendidas durante el curso de Biopython de la licenciatura de ciencias genomicas,cursado durante el tercer semestre, perteneciente a la generacion 2021, el proposito principal fue representar las habilidades aprendidas en el curso, ademas de esto se ha desarrollado un proyecto interno en relacion a cupriavidus, un genero de bacterias que recientemente ha ingerido a otro genero como lo fue Wautersia, por lo que nos interesa saber su similitud, determinando su tamano de genoma, la propia taxonomia; asi como una serie de datos mas como los centros de investigacion donde se han secuenciado las bacterias o tecnologias por las que se han secuenciado.
Un poco sobre cupriavidus:

**En que consiste el genero Cupriavidus?**

Cupriavidus es un genero de bacterias, estas bacterias se caracterizan por ser del tipo Gram-negativas; normalmente tienen forma de “baston” y poseen flagelos, ademas de ser organismos aeróbicos obligatorios, y son quimioorganotróficos o quimiolitotroficos.   El género Cupriavidus fue propuesto en 1987 por Makkar y Casida quienes describieron un depredador bacteriano no obligado de bacterias en el suelo como Cupriavidus necator ( Makkar & Casida, 1987), que fue la primera especie encontrada del genero Cupriavidus, este género está formado por 17 especies distintas, tomando en cuenta que Wautersia (un género que se integró a Cupriavidus en 2005, porque se determinó que era igual al genero Cupriavidus ya existente), asi que no es un genero extremadamente grande. Cabe aclarar que son 17 especies reconocidas por la comunidad cientifica, pero que es muy  complicado diferenciar especies dentro del dominio de las bacterias (diferenciadas principalmente por  las secuencias del gen de ARNr 16S) así que si bien tenemos más estudios que intentan el reconocimiento de nuevas especies, no han sido oficializadas.


**Cuál es el ambiente optimo de este genero bacteriano?** 

Miembros del genero Cupriavidus han sido aislados de diversos ambientes ecologicos como suelo, agua, sedimentos de estanques, nodulos de leguminosas, depositos de flujo de lodo volcanico y fuentes clínicas humanas. 


**Y por que seria de interes este genero para considerar su estudio?** 

Se han encontrado una variedad de funciones fisiologicas en los miembros del genero Cupriavidus; por ejemplo, Cupriavidus taiwanensis puede nodular diferentes especies de leguminosas (da Silva et al, 2012), ademas  algunas especies poseen características de promoción del crecimiento de las plantas, como la secreción de ácido indol acetico y sideroforos y la solubilizacion del fosfato (Pongsilp et al,2012), y aun mas interesante se ha informado que Cupriavidus pampae y C. numazuensis degradan contaminantes organicos como el acido 4- (2,4-diclorofenoxi) butírico (2,4-DB) y el tricloroetileno.

## Metodologia 
Para poder obtener las comparaciones de especies de Cupriavidus, desarrollamos una paqueteria la cual contiene un modulo al que se le dio el nombre de bioanalyser, el cual permite generar principalmente analisis de datos estadisticos obtenidos de ensambles. Este modulo contiene distintas funciones y una unica clase, las cuales obtienen su informacion de NCBI.

En nuestro analisis empleamos unicamente las funciones stats_dataframe y stats_graph.

Con la funcion stats_dataframe generamos un archivo con extension csv el cual comparo las estadisticas de los reportes de los ensambles de todas las especies de Cupriavidus, excluyendo a los Cupriavidus no clasificados y a los provenientes de muestras ambientales.
Mediante la funcion stats_graph generamos graficas de los valores de las estadisticas, donde comparamos el numero de ocurrencias para las estadisticos que no eran valores numericos, y comparamos el valor de los estadisticos que si eran numericos. Los estadisticos que comparamos fueron los siguientes: longitud total del ensamble, nombres de los organismos, nivel del ensamble, representacion del genoma, metodos de ensamble, cobertura del genoma, tecnologias de secuenciacion, y los submitters

Para generar estas graficas requerimos el uso de los siguientes comandos desde terminal
python bioanalyser.py -f 4 -g "Cupriavidus agavae,Cupriavidus alkaliphilus,Cupriavidus basilensis,Cupriavidus campinensis,Cupriavidus cauae,Cupriavidus gilardii,Cupriavidus lacunae,Cupriavidus laharis,Cupriavidus malaysiensis,Cupriavidus metallidurans,Cupriavidus nantongensis,Cupriavidus necator,Cupriavidus neocaledonicus,Cupriavidus numazuensis,Cupriavidus oxalaticus,Cupriavidus pampae,Cupriavidus pauculus,Cupriavidus pinatubonensis,Cupriavidus plantarum,Cupriavidus respiraculi,Cupriavidus taiwanensis,Cupriavidus yeoncheonensis" -e iramirez@lcg.unam.mx -o images/cupriavidus_assemblies_total_length.png -s "total-length" -w 30 -t 10

python bioanalyser.py -f 4 -g "Cupriavidus agavae,Cupriavidus alkaliphilus,Cupriavidus basilensis,Cupriavidus campinensis,Cupriavidus cauae,Cupriavidus gilardii,Cupriavidus lacunae,Cupriavidus laharis,Cupriavidus malaysiensis,Cupriavidus metallidurans,Cupriavidus nantongensis,Cupriavidus necator,Cupriavidus neocaledonicus,Cupriavidus numazuensis,Cupriavidus oxalaticus,Cupriavidus pampae,Cupriavidus pauculus,Cupriavidus pinatubonensis,Cupriavidus plantarum,Cupriavidus respiraculi,Cupriavidus taiwanensis,Cupriavidus yeoncheonensis" -e iramirez@lcg.unam.mx -o images/cupriavidus_assemblies_organism_name.png -s "Organism name" -w 30 -t 10

python bioanalyser.py -f 4 -g "Cupriavidus agavae,Cupriavidus alkaliphilus,Cupriavidus basilensis,Cupriavidus campinensis,Cupriavidus cauae,Cupriavidus gilardii,Cupriavidus lacunae,Cupriavidus laharis,Cupriavidus malaysiensis,Cupriavidus metallidurans,Cupriavidus nantongensis,Cupriavidus necator,Cupriavidus neocaledonicus,Cupriavidus numazuensis,Cupriavidus oxalaticus,Cupriavidus pampae,Cupriavidus pauculus,Cupriavidus pinatubonensis,Cupriavidus plantarum,Cupriavidus respiraculi,Cupriavidus taiwanensis,Cupriavidus yeoncheonensis" -e iramirez@lcg.unam.mx -o images/cupriavidus_assemblies_assembly_level.png -s "Assembly level" -w 30 -t 10

python bioanalyser.py -f 4 -g "Cupriavidus agavae,Cupriavidus alkaliphilus,Cupriavidus basilensis,Cupriavidus campinensis,Cupriavidus cauae,Cupriavidus gilardii,Cupriavidus lacunae,Cupriavidus laharis,Cupriavidus malaysiensis,Cupriavidus metallidurans,Cupriavidus nantongensis,Cupriavidus necator,Cupriavidus neocaledonicus,Cupriavidus numazuensis,Cupriavidus oxalaticus,Cupriavidus pampae,Cupriavidus pauculus,Cupriavidus pinatubonensis,Cupriavidus plantarum,Cupriavidus respiraculi,Cupriavidus taiwanensis,Cupriavidus yeoncheonensis" -e iramirez@lcg.unam.mx -o images/cupriavidus_assemblies_assembly_method.png -s "Assembly method" -w 30 -t 10

python bioanalyser.py -f 4 -g "Cupriavidus agavae,Cupriavidus alkaliphilus,Cupriavidus basilensis,Cupriavidus campinensis,Cupriavidus cauae,Cupriavidus gilardii,Cupriavidus lacunae,Cupriavidus laharis,Cupriavidus malaysiensis,Cupriavidus metallidurans,Cupriavidus nantongensis,Cupriavidus necator,Cupriavidus neocaledonicus,Cupriavidus numazuensis,Cupriavidus oxalaticus,Cupriavidus pampae,Cupriavidus pauculus,Cupriavidus pinatubonensis,Cupriavidus plantarum,Cupriavidus respiraculi,Cupriavidus taiwanensis,Cupriavidus yeoncheonensis" -e iramirez@lcg.unam.mx -o images/cupriavidus_assemblies_genome_coverage.png -s "Genome coverage" -w 30 -t 10

python bioanalyser.py -f 4 -g "Cupriavidus agavae,Cupriavidus alkaliphilus,Cupriavidus basilensis,Cupriavidus campinensis,Cupriavidus cauae,Cupriavidus gilardii,Cupriavidus lacunae,Cupriavidus laharis,Cupriavidus malaysiensis,Cupriavidus metallidurans,Cupriavidus nantongensis,Cupriavidus necator,Cupriavidus neocaledonicus,Cupriavidus numazuensis,Cupriavidus oxalaticus,Cupriavidus pampae,Cupriavidus pauculus,Cupriavidus pinatubonensis,Cupriavidus plantarum,Cupriavidus respiraculi,Cupriavidus taiwanensis,Cupriavidus yeoncheonensis" -e iramirez@lcg.unam.mx -o images/cupriavidus_assemblies_sequencing_technology.png -s "Sequencing technology" -w 30 -t 10

python bioanalyser.py -f 4 -g "Cupriavidus agavae,Cupriavidus alkaliphilus,Cupriavidus basilensis,Cupriavidus campinensis,Cupriavidus cauae,Cupriavidus gilardii,Cupriavidus lacunae,Cupriavidus laharis,Cupriavidus malaysiensis,Cupriavidus metallidurans,Cupriavidus nantongensis,Cupriavidus necator,Cupriavidus neocaledonicus,Cupriavidus numazuensis,Cupriavidus oxalaticus,Cupriavidus pampae,Cupriavidus pauculus,Cupriavidus pinatubonensis,Cupriavidus plantarum,Cupriavidus respiraculi,Cupriavidus taiwanensis,Cupriavidus yeoncheonensis" -e iramirez@lcg.unam.mx -o images/cupriavidus_assemblies_submitter.png -s "Submitter" -w 30 -t 10

Para la generacion del archivo csv que es una tabla de todas las estadisticas, se requirio el siguiente comando
python bioanalyser.py -f 3 -g "Cupriavidus agavae,Cupriavidus alkaliphilus,Cupriavidus basilensis,Cupriavidus campinensis,Cupriavidus cauae,Cupriavidus gilardii,Cupriavidus lacunae,Cupriavidus laharis,Cupriavidus malaysiensis,Cupriavidus metallidurans,Cupriavidus nantongensis,Cupriavidus necator,Cupriavidus neocaledonicus,Cupriavidus numazuensis,Cupriavidus oxalaticus,Cupriavidus pampae,Cupriavidus pauculus,Cupriavidus pinatubonensis,Cupriavidus plantarum,Cupriavidus respiraculi,Cupriavidus taiwanensis,Cupriavidus yeoncheonensis" -e iramirez@lcg.unam.mx -o data/Cupriavidus_Assemblies_Stats.csv

## Pregunta 


**Podemos obtener información de las distintas caracteristicas entre las diversas especies de Cupriavidus mediante el uso de librerias de Python?**
El objetivo a niveles generales es obtener la mayor cantidad posible de informacion que nos puede ofrecer NCBI sobre el genero cupriavidus. 


## Resultados  
Dejando en claro que para la informacion presente en estos resultados no se consideraron las Cupriavidus no clasificadas y las de muestra ambiente; asi que no influyen en la informacio presentada dentro de este contenido. 
Los datos de los resultados principales pueden observarse de manera ordenada dentro de nuestro apartado de Data en un archivo CSV, ahi se contiene informacion esencial como lo fue el numero total de nucleotidos del genoma, asi como si este se encuentra cerrado o abierto, el porcentaje de CG, tambien el numero de gaps y scaffolds, asi como la tecnologia por la cual fue secuenciada, el ano de publicacion o el instituto o centro de investigacion que ha publicado la informacion; claro tambien una de los aspectos mas importantes como lo son las especies especificas dentro del genero.  

**Visualizando el tamano de genoma para cada uno de los organismos secuenciados e identificados**

Dentro de los resultados hemos observado por ejemplo una gran diferencia entre el genoma mas pequeno del genero y el mas grande, siendo estos: ASM782943v1 como el mas pequeno con un total de 4879332 nucleotidos siendo esta perteneciente a la especie Cupriavidus gilardii y el mas grande ha sido ASM74409v1 perteneciente a la especie Cupriavidus necator de la cepa A5-1; como podemos observar es una diferencia gigantesca de 4862044 nucleotidos o sea que es practicamente el doble de grande, esto es muy llamativo ya que como ya hemos mencionado las bacterias son apenas diferenciables dentro del genero, pero en este caso son bastante distintas.
![image](https://user-images.githubusercontent.com/60999318/144784947-92d1fbb6-550d-4651-b9ef-17f2274bad66.png)


**Numero de representaciones de especies cupriavidus**

Las especies con mayor numero de representantes identificados han sido taiwanensis, metalladurans y necator, esto debe ser por su importancia en ambitos de investigacion reciente como puede ser la biorremediacion, o como especies simbioticas con plantas. 

![image](https://user-images.githubusercontent.com/60999318/144790299-7236ec7c-4023-4f23-8e55-f39184220532.png)


**Representacion del nivel de ensamblaje**

Como podemos observar en esta imagen la distribucion del nivel de ensamblaje es muy poco distante, sin embargo el mayor es de tipo contig, aunque cromosomas y genomas completos estan muy presentes.

![image](https://user-images.githubusercontent.com/60999318/144790820-c33d5a94-c203-448d-8af5-16da7765ce3f.png)

**Institutos que contribuyeron a la secuenciacion de Cupriavidus**

Dentro del top 5 de este aspecto tenemos al instituto Pasteur, DOE joint Genome Institute, LM-UGent, US food and drug administration y el CCUG 
![image](https://user-images.githubusercontent.com/60999318/144791127-5c0973ad-2f2b-4500-9683-e4acabffa9e6.png)

 **Metodos de ensamblaje**
 
 nan que es el visiblemente mas numero, realmente lo que nos dice es que los investigadores no proporcionaron esta informacion, sin embargo de la informacion que si conocemos es que Spades es el ensamblador mas comun a la hora de ensamblar los genomas de cupriavidus. 
 
 ![image](https://user-images.githubusercontent.com/60999318/144791521-de433add-6e82-437b-a633-8c898e0d4fcf.png)


**Tecnologia de secuenciacion usada**

Aqui podemos observar que nuevamente, por lo regular no se reporta la tecnologia de secuenciacion usada para la obtencion de estos genomas (algo bastante malo), sin embargo, dentro de la informacion que si nos es proporcionada podemos determinar que Illumina y sus variantes es la tecnologia preferida a la hora de secuenciar este genero. 

![image](https://user-images.githubusercontent.com/60999318/144792026-0905d611-64fe-4299-a785-b74d85a48628.png)


**Cobertura de genoma**
Observamos que existe una gran variacion de la precision de la secuenciacion entre los distintos organismos de Cupriavidus.

![image](https://user-images.githubusercontent.com/60999318/144794891-d6cc7f9b-b73f-4575-9919-498d4d12c888.png)



## Conclusiones
El uso de herramientas bioinformaticas ha sido lo que nos ha permitido dar un salto gigantesco en el desarrollo de la ciencia, en este caso hablamos en el ambito de la biologia, pero sin duda alguna la ciencia de datos ha tenido un impacto gigante en toda la ciencia, desde la fisica hasta la astronomia, por lo que es importante que aprendamos a interpretar el mundo que nos rodea a traves del computo; quiza un sin fin  de numeros, letras y codigo no nos pueda decir demasiado, pero somos capaces de sintetizar todo este mar de informacion dentro de graficas o tablas.  Tambien reconocer que el buen manejo de estas herramientas, la capacidad de acceder a la informacion e interpretarla nos volvera unos mejores cientificos a futuro, personalmente tambien el desarrollo de este proyecto me ha ensanado la importancia del orden y la anticipacion.  
 En el apartado de Cupriavidus podemos concluir que es un genero que aunque apenas esta comenzando a poseer relevancia es muy diverso e interesante, casos como lo son la diferencia entre el tamano de genoma de su miembro mayor y menor nos han ensenado que la familiaridad de estos organismos depende de mucho mas que tan solo de un genoma casi identico. Sin duda aun tenemos mucho por descubrir de este genero tan poco estudiado. 
 Hoy en dia tenemos la venteja de que la informacion esta al alcance de nuestras manos, por lo que seria un insulto no ser capaces de aprovechar esta ventaja y este contexto historico para que nosotros tambien nos volvamos un cumulo de informacion.
