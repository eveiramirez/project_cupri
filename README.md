# project_cupri


## INTRODUCCION
En este proyecto se trabajo con las herramientas aprendidas durante el curso de Biopython de la licenciatura de ciencias genomicas,cursado durante el tercer semestre, perteneciente a la generacion 2021, el proposito principal fue representar las habilidades aprendidas en el curso, ademas de esto se ha desarrollado un proyecto interno en relacion a cupriavidus, un genero de bacterias que recientemente ha ingerido a otro genero como lo fue Wautersia, por lo que nos interesa saber su similitud, determinando su tamano de genoma, la propia taxonomia; asi como una serie de datos mas como los centros de investigacion donde se han secuenciado las bacterias o tecnologias por las que se han secuenciado.
Un poco sobre cupriavidus:

**En qué consiste el género Cupriavidus?**

Cupriavidus es un genero de bacterias, estas bacterias se caracterizan por ser del tipo Gram-negativas; normalmente tienen forma de “baston” y poseen flagelos, ademas de ser organismos aeróbicos obligatorios, y son quimioorganotróficos o quimiolitotroficos.   El género Cupriavidus fue propuesto en 1987 por Makkar y Casida quienes describieron un depredador bacteriano no obligado de bacterias en el suelo como Cupriavidus necator ( Makkar & Casida, 1987), que fue la primera especie encontrada del genero Cupriavidus, este género está formado por 17 especies distintas, tomando en cuenta que Wautersia (un género que se integró a Cupriavidus en 2005, porque se determinó que era igual al genero Cupriavidus ya existente), asi que no es un genero extremadamente grande. Cabe aclarar que son 17 especies reconocidas por la comunidad cientifica, pero que es muy  complicado diferenciar especies dentro del dominio de las bacterias (diferenciadas principalmente por  las secuencias del gen de ARNr 16S) así que si bien tenemos más estudios que intentan el reconocimiento de nuevas especies, no han sido oficializadas.


**Cuál es el ambiente óptimo de este género bacteriano?** 

Miembros del genero Cupriavidus han sido aislados de diversos ambientes ecologicos como suelo, agua, sedimentos de estanques, nodulos de leguminosas, depositos de flujo de lodo volcanico y fuentes clínicas humanas. 


**Y por qué sería de interés este género para considerar su estudio?** 

Se han encontrado una variedad de funciones fisiologicas en los miembros del genero Cupriavidus; por ejemplo, Cupriavidus taiwanensis puede nodular diferentes especies de leguminosas (da Silva et al, 2012), ademas  algunas especies poseen características de promoción del crecimiento de las plantas, como la secreción de ácido indol acetico y sideroforos y la solubilizacion del fosfato (Pongsilp et al,2012), y aun mas interesante se ha informado que Cupriavidus pampae y C. numazuensis degradan contaminantes organicos como el acido 4- (2,4-diclorofenoxi) butírico (2,4-DB) y el tricloroetileno.

## Metodologia


## Pregunta


**Podemos obtener información de las distintas características entre las diversas especies de Cupriavidus mediante el uso de librerías de Python?**
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
