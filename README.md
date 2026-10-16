# Proyecto Final - Transcriptómica 
### Héctor Jesús López-Ordaz
#### 2026-05-31

___

# Planteamiento
La fibrosis pulmonar idiopática (IPF) es una enfermedad pulmonar progresiva de etiología desconocida, caracterizada por una fibrosis pulmonar progresiva [1]. Caracterizar los cambios transcriptómicos entre el estado sano y de IPF ofrece la posibilidad de entender a detalle qué procesos moleculares son responsables del cambio de la arquitectura pulmonar. 
Con el fin de aplicar los conocimientos adquiridos, se realizó el siguiente proyecto, que busca cubrir dos objetivos; en primera instancia, aplicar los conocimientos adquiridos a lo largo del curso de Transcriptómica en combinación al entendimiento de los procesos biológicos, esto con el fin de comprender que la bioinformática y la biología se ven enriquecidas cuando las disciplinas son integradas. En segunda instancia comprender las diferencias que caracterizan a la IPF contra un control sano, buscando que el conocimiento adquirido concuerde con literatura y comprensión propia del tema.
Para llevar a cabo los puntos antes descritos, se decidió usar datos provenientes de un artículo de enfoque multiómico que tuvo el objetivo de identificar las señales iniciales de la IPF [1], disponible en Gene Expression Omnibus (GEO) con número de acceso [GSE213001](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213001). El artículo proveía de los datos específicos que se buscó indagar (muestras de IPF y sanas), contenía sección de reproducibilidad que facilitaba tratar los datos de manera fiel, mostraba gráficos de componentes principales que anticipa la de diferenciación entre grupos de análisis, y como extra el número de acceso de GEO agrupa una cantidad masiva de archivos de secuencias de lectura (SRR) lo cuál consolidó la elección. Del set total se seleccionó un subset con un filtro. El filtro consistió para el caso de IPF, muestras de la base del pulmón excluyendo el ápice pulmonar, y para el caso del grupo de control sin enfermedad (NDC), no hubo discriminación por localización, esto debido a que el mismo estudio ofrece pruebas estadísticas que demuestran la no diferencia entre localización de las muestras NDC. Además para evitar la pseudoreplicación, las muestras de IPF y NDC provienen de distintos pacientes. El filtro arrojó catorce muestras NDC y veinte IPF.
Las herramientas utilizadas para el análisis fueron aquellas creadas alrededor del curso y que se describe en cada reporte enviado como tarea. De manera breve siguen un flujo bastante directo donde la entrada de uno es la salida del anterior. La colección de herramientas incluye el uso de programas para análisis de control de calidad de SRR (`fastqc`, `multiqc`), limpieza de las mismas (`fastp`), alineamiento contra el genoma de referencia (`STAR`), construcción de *feature counts* (`featureCounts`), análisis de expresión diferencial (`DESeq2`), anotación funcional (DAVID y GSEA) e interacción de proteína-proteína (IPP) (STRING). Todo esto con el uso crítico de las *core utils* de Linux y R para manejar el flujo, organización e incluso creación de los archivos generados.
# Procedimiento y Justificación

> Todos los datos se encuentran en el servidor “chaac” de la Licenciatura bajo la ruta:
> `/home/hectorjl/4to/transcriptomica_proyecto_final`

Como primer paso se accedió al GSE213001 desde GEO, donde se capturó el ID del Bioproject y se procedió a buscar en (SRA Run Selector)[https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA878759&o=acc_s%3Aa], desde ahí se descargó la tabla de metadatos (`SraRunTable.csv`).

## Selección de SRRs

Se utilizó las *core utils* para filtrar los SRRs que se descargaron (comando en Anexo n.1). El filtro seleccionaba muestras IPF de la base ó NDC de la base o ápice del pulmón, además que para ambos casos el paciente no debería repetirse, asegurando muestras independientes. Catorce muestras NDC y veinte IPF pasaron el filtro y se muestran a continuación:

| Run         | donorid | diagnosis |
|:-:|:-:|:-:|
| SRR21498285 | ALF037  | IPF       |
| SRR21498295 | ALF032  | IPF       |
| SRR21498375 | ALF019  | IPF       |
| SRR21498185 | ALF027  | IPF       |
| SRR21498421 | ALF028  | IPF       |
| SRR21498828 | ALF015  | IPF       |
| SRR21498668 | ALF035  | IPF       |
| SRR21498365 | ALF029  | IPF       |
| SRR21498165 | ALF023  | IPF       |
| SRR21498205 | ALF001  | IPF       |
| SRR21498628 | ALF018  | IPF       |
| SRR21498280 | ALF045  | IPF       |
| SRR21498230 | ALF008  | IPF       |
| SRR21498643 | ALF041  | IPF       |
| SRR21498200 | ALF039  | IPF       |
| SRR21498748 | ALF016  | IPF       |
| SRR21498340 | ALF043  | IPF       |
| SRR21498265 | ALF047  | IPF       |
| SRR21498538 | ALF040  | IPF       |
| SRR21498633 | ALF038  | IPF       |
| SRR21498578 | ALF013  | NDC       |
| SRR21498190 | ALF026  | NDC       |
| SRR21498195 | ALF030  | NDC       |
| SRR21498663 | ALF009  | NDC       |
| SRR21498225 | ALF020  | NDC       |
| SRR21498150 | ALF050  | NDC       |
| SRR21498523 | ALF024  | NDC       |
| SRR21498320 | ALF025  | NDC       |
| SRR21498623 | ALF012  | NDC       |
| SRR21498400 | ALF044  | NDC       |
| SRR21498140 | ALF049  | NDC       |
| SRR21498385 | ALF034  | NDC       |
| SRR21498728 | ALF010  | NDC       |
| SRR21498180 | ALF017  | NDC       |

**Tabla 1.** Muestras seleccionadas para el análisis con su correspondiente SRR, se incluye también el ID del donante que no se repite en ningún caso.

## Descarga de SRRs y conversión a archivos FASTQ

```bash
# fasterq-dump:
#	--outdir: Directorio de salida donde los fastq descargados se guardan
prefetch (cat data/GSE213001/SRR_list.txt) && for line in $(cat data/GSE213001/SRR_list.txt); do fasterq-dump --outdir data/raw_fastqs "$line" done
```

Tomando los SRRs elegidos, se procedió a la descarga de los archivos intermedios SRA y con ellos los archivos finales FASTQ. La combinación de `prefetch` y `fasterq-dump` es la forma más rápida de extraer archivos FASTQ de los registros de SRA [2].
**Resultados**: La descarga finalizó con treinta y cuatro archivos FASTQ debido a que las lecturas fueron de tipo *single end* correspondientes a catorce de muestras NDC y veinte de IPF. Los archivos FASTQ fueron renombrados usando el comando del Anexo n.2

## Control de calidad


```bash
fastqc -t 34 data/raw_fastqs/*.fastq -o data/raw_qcs 
multiqc data/raw_qcs/* -o data/raw_multiqc
```

FastQC es una aplicación que toma un archivo FASTQ y le aplica una serie de pruebas para generar un informe completo de control de calidad. Este informe indica si hay algo inusual en las lecturas [3]. Mientras MultiQC crea un compendio global de todos los reportes generados por FastQC. Los informes generados resultan fundamentales para la creación de umbrales de limpieza para el posterior procesamiento.
**Resultados**:FastQC generó un informe HTML para cada lectura. Se revisó uno a uno y se concluyó ciertos aspectos que aplicaban para todas las lecturas. Se usó el reporte integrador HTML de MultiQC para reafirmar los hallazgos individuales. Con la información adquirida de las herramientas se concluyó lo siguiente: de manera global la calidad de la secuenciación fue muy alta (calidad por base > 32), con una longitud de cien bases por lectura; las únicas métricas que fueron de notar incluyeron la presencia de adaptadores (poly a e Illumina) y alteración en la proporción de nucleótidos de las doce primeras bases de cada lectura.

## Limpieza de lecturas y control de calidad post-limpieza

```bash
# Código fuente completo: src/auto_fastp.sh
# -- [recorte] --
# 		Parámetros:
   fastp -i "$1" \ # Path de archivo FASTQ de entrada
          -o "$2" \ # Path de archivo FASTQ de salida
          --trim_front1 12 \ #  Recorta las primeras 12 bases de cada lectura
          --trim_poly_g \ # Recorta colas poly g, introducidas por errores del secuenciador
          --trim_poly_x \ # Recorta colas x (a, c, g), introducidas por error del secuenciador
-l 50 # Elimina lecturas si no cumplen con una longitud mínima de 50 bases
# -- [recorte] --
```
Para la limpieza se usó `fastp`, que es una herramienta ultrarrápida para realizar control de calidad, el filtrado de lecturas y la corrección de bases en datos FASTQ [4]. La elección de los argumentos usados en su ejecución fue basada en las observaciones del paso previo, buscando reducir el ruido que podrían introducir más adelante en ejecuciones posteriores.
Se ejecutó sobre las lecturas procesadas nuevamente el control de calidad con FastQC y MultiQC con el fin de comprobar que el procesamiento realizó lo previsto.
**Resultados**: Se generó un total de treinta y cuatro lecturas “limpias”. Los reportes de QC mostraron una mejor calidad global, adaptadores de Illumina  fueron eliminados y de poly a reducidos, además la alteración en la proporción de nucleótidos de las primeras bases fué eliminada. Se presenta la comparación global (lecturas crudas contra procesadas) de estadísticas generadas por MultiQC:



![Heatmap de métricas crudas](https://raw.githubusercontent.com/Hectorjlo/transcriptomica_proyecto_final/refs/heads/main/data/raw_multiqc/plots/raw_fastqc-status-check-heatmap.png)


![Heatmap de métricas limpias](https://raw.githubusercontent.com/Hectorjlo/transcriptomica_proyecto_final/refs/heads/main/data/cleaned_multiqc/plots/cleaned_fastqc-status-check-heatmap.png)



Se concluyó que la limpieza fue exitosa y posibilitó continuar para el próximo paso.

## Alineamiento de lecturas contra genoma de referencia.
Antes de la alineación, se realizó la creación del índice necesario para el alineador, este fué creado con el genoma de referencia humano [GRCh38.104](https://www.gencodegenes.org/human/releases.html) (comando en Anexo n.3) obtenido de GENCODE, en concordancia al usado en el estudio donde se obtuvo los datos (artículo de referencia). Por último los archivos SAM fueron ordenados y transformados en BAM.

```bash
# -- [recorte] --
# Código fuente completo: src/auto_star_single.sh
# 	# Parámetros:
   STAR --runMode alignReads \ # Modo de alinear lecturas contra un índice
	# Path al índice del genoma
         --genomeDir /home/hectorjl/WORKING_DIR/4to/transcriptomica_proyecto_final/data/GENCODE_GRCh38.p13_104/STAR_index \
         --readFilesIn "$1" \ # Path al archivo de lecturas a alinear
         --outSAMtype SAM \ # Formato de salida del alineamiento
         --outFileNamePrefix "$2" \ # Path y prefijo de los archivos generados por STAR
         --outSAMunmapped None \ # Elimina entradas de genes no mapeados por la lectura
         --runThreadN 8 \ # Utiliza 8 hilos por proceso
         --genomeLoad LoadAndKeep \ # Carga el índice desde un espacio de memoria evitando cargar un genoma por proceso
         --outFilterMultimapNmax 7 # Cantidad máxima de alineamientos al genoma por lectura
# -- [recorte] --
```
El uso de STAR (Spliced Transcripts Alignment to a Reference) se debe en parte a que fue el usado en el artículo de referencia y a que es un alineador universal ultrarrápido para secuenciación de RNA [5]. Además que a lo largo del curso fue la herramienta que personalmente proporcionó mejor alineación global.

**Resultados**: En concordancia con la calidad de las lecturas, STAR pudo alinear el más del 78% de las lecturas al genoma (Mapeo único). Se resumen los resultados a continuación:

| Muestra (SRR) | Grupo | Mapped Multiple Loci | Uniquely Mapped Reads | Suma |
|:-:|:-:|:-:|:-:|:-:|
| SRR21498165          | IPF   | 11.45%                 | 86.64%                  | 98.09%   |
| SRR21498185          | IPF   | 10.64%                 | 87.69%                  | 98.33%   |
| SRR21498200          | IPF   | 10.81%                 | 86.95%                  | 97.76%   |
| SRR21498205          | IPF   | 16.58%                 | 80.84%                  | 97.42%   |
| SRR21498230          | IPF   | 11.63%                 | 85.89%                  | 97.52%   |
| SRR21498265          | IPF   | 17.96%                 | 78.45%                  | 96.41%   |
| SRR21498280          | IPF   | 15.89%                 | 81.48%                  | 97.37%   |
| SRR21498285          | IPF   | 10.25%                 | 88.04%                  | 98.29%   |
| SRR21498295          | IPF   | 11.19%                 | 86.63%                  | 97.82%   |
| SRR21498340          | IPF   | 9.58%                  | 88.82%                  | 98.40%   |
| SRR21498365          | IPF   | 17.44%                 | 79.90%                  | 97.34%   |
| SRR21498375          | IPF   | 14.46%                 | 83.18%                  | 97.64%   |
| SRR21498421          | IPF   | 10.84%                 | 87.49%                  | 98.33%   |
| SRR21498538          | IPF   | 12.56%                 | 85.05%                  | 97.61%   |
| SRR21498628          | IPF   | 11.59%                 | 86.44%                  | 98.03%   |
| SRR21498633          | IPF   | 10.09%                 | 88.41%                  | 98.50%   |
| SRR21498643          | IPF   | 9.64%                  | 88.79%                  | 98.43%   |
| SRR21498668          | IPF   | 10.00%                 | 88.29%                  | 98.29%   |
| SRR21498748          | IPF   | 10.31%                 | 88.07%                  | 98.38%   |
| SRR21498828          | IPF   | 11.58%                 | 86.10%                  | 97.68%   |
| SRR21498140          | NDC   | 10.47%                 | 86.66%                  | 97.13%   |
| SRR21498150          | NDC   | 9.69%                  | 88.28%                  | 97.97%   |
| SRR21498180          | NDC   | 11.20%                 | 86.79%                  | 97.99%   |
| SRR21498190          | NDC   | 9.96%                  | 88.31%                  | 98.27%   |
| SRR21498195          | NDC   | 8.94%                  | 89.22%                  | 98.16%   |
| SRR21498225          | NDC   | 9.53%                  | 88.39%                  | 97.92%   |
| SRR21498320          | NDC   | 10.36%                 | 87.88%                  | 98.24%   |
| SRR21498385          | NDC   | 8.98%                  | 89.68%                  | 98.66%   |
| SRR21498400          | NDC   | 10.02%                 | 87.85%                  | 97.87%   |
| SRR21498523          | NDC   | 9.99%                  | 88.04%                  | 98.03%   |
| SRR21498578          | NDC   | 10.31%                 | 87.12%                  | 97.43%   |
| SRR21498623          | NDC   | 10.83%                 | 87.34%                  | 98.17%   |
| SRR21498663          | NDC   | 10.01%                 | 87.88%                  | 97.89%   |
| SRR21498728          | NDC   | 10.02%                 | 87.65%                  | 97.67%   |

**Tabla 2:** Porcentajes de mapeo de las lecturas contra el genoma de referencia.

Se obtuvo un archivo SAM por cada alineamiento, que fue transformado en archivos tipo BAM ordenados y más ligeros (comando en Anexo n.4). Los altos porcentajes de mapeo indican que las lecturas corresponden correctamente al genoma de referencia usado. Lo que ofrece mejores resultados y da más rigurosidad estadística en los próximos pasos

## Mapeo de cuentas por gen y creación de matriz de conteos

```bash
# -- [recorte] –
# Código fuente completo: src/auto_ft_counts.sh
#			Parámetros:
featureCounts -o "$OUTPUT_FILE" \ # Path al archivo de salida
    -a "$ANNOTATION" \ # Archivo GFF3 de anotación del genoma de referencia
    -T 24 \ # Cantidad de hilos para usar durante el proceso
	# Si hay superposición de una lectura a distintos features, asigna la lectura al feature con el que comparte más bases
    --largestOverlap \ 
    "${bam_files[@]}" # Array de Paths a los archivos de entrada
# -- [recorte] --
```

Para hacer el análisis de expresión diferencial primero se necesita obtener las cuentas por gen. Para ello se utiliza `featureCounts`, programa de resumen de lecturas de gran eficiencia que cuenta las lecturas mapeadas a características genómicas como exones [6]. Se utilizó su salida para crear la matriz de conteos y el archivo que contiene el ID del gene y su longitud, que son necesarios para el análisis de expresión diferencial.
**Resultados**: Del stdout de la ejecución se recuperó la cantidad porcentual de lecturas asignadas a features, y se resumen en la siguiente tabla:
| Muestra (SRR)       | Grupo | Alineamientos Asignados |
|:-:|:-:|:-:|
| SRR21498140          | NDC   | 90.5%                     |
| SRR21498150          | NDC   | 87.4%                     |
| SRR21498180          | NDC   | 88.5%                     |
| SRR21498190          | NDC   | 89.2%                     |
| SRR21498195          | NDC   | 86.2%                     |
| SRR21498225          | NDC   | 88.0%                     |
| SRR21498320          | NDC   | 88.5%                     |
| SRR21498385          | NDC   | 88.2%                     |
| SRR21498400          | NDC   | 87.2%                     |
| SRR21498523          | NDC   | 87.9%                     |
| SRR21498578          | NDC   | 89.0%                     |
| SRR21498623          | NDC   | 91.8%                     |
| SRR21498663          | NDC   | 91.7%                     |
| SRR21498728          | NDC   | 91.3%                     |
| SRR21498165          | IPF   | 86.7%                     |
| SRR21498185          | IPF   | 89.7%                     |
| SRR21498200          | IPF   | 83.4%                     |
| SRR21498205          | IPF   | 91.6%                     |
| SRR21498230          | IPF   | 88.1%                     |
| SRR21498265          | IPF   | 86.3%                     |
| SRR21498280          | IPF   | 87.3%                     |
| SRR21498285          | IPF   | 85.3%                     |
| SRR21498295          | IPF   | 86.2%                     |
| SRR21498340          | IPF   | 86.5%                     |
| SRR21498365          | IPF   | 90.3%                     |
| SRR21498375          | IPF   | 88.3%                     |
| SRR21498421          | IPF   | 87.7%                     |
| SRR21498538          | IPF   | 84.9%                     |
| SRR21498628          | IPF   | 85.4%                     |
| SRR21498633          | IPF   | 87.6%                     |
| SRR21498643          | IPF   | 90.5%                     |
| SRR21498668          | IPF   | 87.3%                     |
| SRR21498748          | IPF   | 87.8%                     |
| SRR21498828          | IPF   | 86.9%                     |

**Tabla 3:** Porcentajes de lecturas mapeadas a features por `featureCounts`.

Se observa de nueva cuenta que la cantidad porcentual global es alta, lo que indica que los datos lograron agruparse en transcritos y resultan útiles para el análisis de expresión diferencial. La tabla generada por la ejecución fue usada para la creación de la matriz de conteos y el archivo que contiene gene ID y su longitud (comando en Anexo n.5).
Al finalizar esta sección se contó con los dos de tres archivos necesarios para el análisis de expresión diferencial.


## Análisis de expresión diferencial

Para el análisis de expresión diferencial se utilizó una combinación de herramientas provenientes de diversas paqueterías de R, el principal fue DESeq2, pero también se usó de manera puntual edgeR, además librerías para graficación como ggplot2 y ComplexHeatmap. A continuación se muestran los bloques críticos para el análisis.
Antes de realizar el análisis se creó el último archivo faltante, este contiene el gene ID y el nombre del gene desde el GTF del genoma de referencia (comando en Anexo n.6).

```R
# -- [recorte] –
# Código fuente completo: src/DEA_DESeq2.R
# Fórmula usada a lo largo del código
## “0 +” Elimina el intercepto, por lo cual cada nivel de 
## “condition” y “sex” se codifique como coeficiente independiente
## “age” Es una covariable numérica continua y no un coeficiente
## por lo cual controla las diferencias de edad entre las muestras

# Creación del objeto de DESeq2 desde la matriz de conteos
dds <- DESeqDataSetFromMatrix(
    countData = round(matrix_counts), # Utiliza la matriz de conteos
    colData = meta_data, # Metadatos que ilustran cómo está organizados los conteos
    design = formula # Fórmula para modelar la expresión
)
# -- [recorte] --
## Con el uso de la función de edgeR “filterByExpr()”
## se eliminó aquellos genes que no tenían suficientes conteos
## y que su mantenimiento podría afectar al ajuste por el FDR
keep <- filterByExpr(dds, design)
# Selecciona solo los genes que pasan el filtro
dds <- dds[keep, ]
# -- [recorte] --
# Creación manual de los contrastes
## Se sustrae la diferencia que existe de IPF menos NDC
## De esta manera las diferencias observadas en los próximos gráficos 
## son respecto al NDC, que actua como nivel basal de expresión.
contrasts <- makeContrasts(IPF_vs_NDC = conditionIPF - conditionNDC,
    levels = design
)

# Realiza el análisis de expresión diferencial
## Se usa el contraste anterior
res <- results(dds, contrast = contrasts[ ,"IPF_vs_NDC"])
# -- [recorte] --
```
El uso de DESeq2, se debe a que es un método para el análisis diferencial, que utiliza la estimación por contracción para las dispersiones y los cambios de expresión con el fin de mejorar la estabilidad y la interpretabilidad de las estimaciones [7]. Además que mostró en la tarea correspondiente, ser el algoritmo más robusto y que mejor maneja los cambios de la expresión génica, sin la necesidad de filtros robustos antes de su uso directo. Por otro lado, el artículo de referencia utilizó a “Degust” para su análisis de expresión diferencial, y a pesar de que no realizan directamente el mismo análisis aquí presente, si ofrece la posibilidad de comprobar la reproducibilidad de un enfoque Web como Degust contra uno más tradicional como DESeq2.

**Resultados**: Los archivos que fueron generados por el análisis son diversos e incluyen: Tablas de genes diferencialmente expresados (DEGs), una tabla de *Transcripts per Million*, la cantidad total de DEGs, el objeto de resultado del análisis RDS; y gráficos como: Heatmaps, Gráficos de volcán y PCAs. Por lo cual estos últimos se presenta en las secciones siguientes de manera ordenada junto a su interpretación.

### Análisis visual
#### PCA
Se creó tres gráficos de análisis de componentes principales para observar la agrupación de la condición, sexo y edad.

![PCA por condición](https://raw.githubusercontent.com/Hectorjlo/transcriptomica_proyecto_final/refs/heads/main/results/star/DESeq2/plots/PCA_condition.png)

![PCA por sexo](https://raw.githubusercontent.com/Hectorjlo/transcriptomica_proyecto_final/refs/heads/main/results/star/DESeq2/plots/PCA_sex.png)

![PCA por edad](https://raw.githubusercontent.com/Hectorjlo/transcriptomica_proyecto_final/refs/heads/main/results/star/DESeq2/plots/PCA_age.png)

Como se puede observar los PCAs de sexo y edad no muestran una agrupación clara por estas características, indicando que la varianza entre las muestras no es explicada plenamente por ellas. Por otro lado, el PCA por condición muestra dos clusters diferenciados entre los grupos IPF y NDC, indicando que esta característica podría explicar las diferencias entre las muestras, además entre los mismos grupos las muestras IPF tienen menor cohesión entre sí, lo que puede indicar que los estados IPF son heterogéneos y no necesariamente comparten todas las características entre sí; a diferencia de los NDC que mantienen una mayor cohesión en su mismo grupo.

#### Gráfico de volcán
Para observar los DEGs considerando la relación entre magnitud de cambio (*log2 Fold Change*) en el eje x, y la significancia estadística (*-log10 FDR*) en el eje y, se realizó un gráfico de volcán.

![Vulcano Plot de DEGs. Los 30 DEGs más expresados se encuentran etiquetados. FDR < 0.001, FC = 0.5](https://raw.githubusercontent.com/Hectorjlo/transcriptomica_proyecto_final/refs/heads/main/results/star/DESeq2/plots/vulcano_plot.png)

Visualizando el gráfico se puede empezar a comparar directamente IPF contra NDC, dado que la base es NDC, los DEGs son específicos de las muestras de IPF. 

**DEGs upregulated:** Se observó la presencia de las familias de las metaloproteasas (MMP7, MMP13, ADAMTS14, SPP1) que están implicadas en la remodelación de la matriz extracelular (MEC) y reclutamiento de fibroblastos [8], además existe presencia de la familia de colágenos (COL17A1, COL7A1, COL14A1) indicando una deposición descontrolada de estas fibras que son marca característica de la IPF [9]. También se localizan genes de señalización (LTPB1, CRABP2) implicados en la vía de de TGF-$\beta$ y reparación tisular disfuncional, marcas claras de un proceso fibrótico.

**DEGs downregulated:** Se encontró a RTKN2 fuertemente suprimido, este gen es marca específica de células epiteliales alveolares de tipo 1 (AEC1s), las cuales son las encargadas del intercambio gaseoso dentro de la cavidad alveolar. Que se encuentre como *downregulated* indica que la función pulmonar está disminuida en IPF. Además de SLCO1A2 y ITLN2 que son potenciales biomarcadores de IPF [8].


#### Heatmaps de DEGs
Recordando que DESeq2 se caracteriza por una contracción de la varianza, priorizando genes con cambios consistentes. A continuación se muestra los heatmaps obtenidos del análisis de expresión diferencial.

![Heatmap DEGs](https://raw.githubusercontent.com/Hectorjlo/transcriptomica_proyecto_final/refs/heads/main/results/star/DESeq2/plots/heatmap_all.png)

![Heatmap top 30 DEGs](https://raw.githubusercontent.com/Hectorjlo/transcriptomica_proyecto_final/refs/heads/main/results/star/DESeq2/plots/heatmap_top30.png)

Recapitulando la heterogeneidad encontrada entre el grupo de IPF, los heatmaps confirman este comportamiento, el grupo de NDC muestra un parecido más amplio entre el mismo grupo, mientras que el grupo IPF tiene una mayor variabilidad visible esto es más presente en el heatmap que incluye todos los DEGs.
Por otro lado, el factor H del complemento (CFH) se muestra *up regulated* en el grupo IPF y que fue el que contó con la mayor significancia estadística de todo el análisis, está implicado en el sistema inmune innato y es un regulador crítico de la patogénesis de la IPF, contribuyendo a la progresión fibrótica marcada de esta enfermedad [10]. 

Por otro lado, se encontró interesante que ARRB1 estuviera *down regulated* en el grupo IPF, esto debido a que la literatura en modelos murinos indica que ARRB1 promueve fibrosis y su deleción protege de la IPF en ratones [11]. Esta diferencia podría explicarse debido a las diferencias entre especies o algún mecanismo aún no descrito, por lo cual abre una posibilidad de investigación en este tema.

## Anotación funcional

Para la anotación funcional, se decidió utilizar la herramienta web [DAVID](https://davidbioinformatics.nih.gov/summary.jsp) (Database for Annotation, Visualization and Integrated Discovery) [12,13]. 
En esta se utilizó como modelo de *background* a los genes que pasaron el filtro (filterByExpr()) del análisis con DESeq2, y como interrogadores a los DEGs *up regulated* y *down regulated*. Después se descargaron los datos los archivos CSV de ontología génica (GO) para *Biological process* (BP) y *Cellular component* (CC) y con ellos los gráficos que se presentan a continuación.

**Resultados:** Se realizó un gráfico para cada set de DEGs y para cada ontología descargada, lo cual acumula a cuatro gráficos totales.
![Enriquecimiento GO downregulated BP](https://raw.githubusercontent.com/Hectorjlo/transcriptomica_proyecto_final/refs/heads/main/results/star/DAVID/plots/GO_Bubble_Plot_downregulated_BP.png)

![Enriquecimiento GO upregulated BP](https://raw.githubusercontent.com/Hectorjlo/transcriptomica_proyecto_final/refs/heads/main/results/star/DAVID/plots/GO_Bubble_Plot_upregulated_BP.png)

![Enriquecimiento GO downregulated CC](https://raw.githubusercontent.com/Hectorjlo/transcriptomica_proyecto_final/refs/heads/main/results/star/DAVID/plots/GO_Bubble_Plot_downregulated_CC.png)

![Enriquecimiento GO upregulated CC](https://raw.githubusercontent.com/Hectorjlo/transcriptomica_proyecto_final/refs/heads/main/results/star/DAVID/plots/GO_Bubble_Plot_upregulated_CC.png)

Los gráficos de BP muestran concordancia con los previos hallazgos reportados en el presente documento; por ejemplo para los BP *down regulated* se muestra un enriquecimiento en la transducción de señales coincidente con el estado celular, donde la MEC rígida afecta a esta transducción. En enriquecimiento en la angiogénesis revela como la fibrosis afecta a la manutención al endotelio alveolar, además como el desorden tisular generado por la fibrosis afecta a la respuesta del sistema inmune como defensa contra patógenos externos. 
Por otro lado BP *up regulated* muestra una historia de enriquecimiento en torno a mayor adhesión celular, producto de la MEC rígida y aumento de integrinas, una respuesta inmune resultado de una inflamación crónica del tejido, y la justa reorganización de la MEC. Los procesos biológicos enriquecidos de cilios y estructuras similares tiene su explicación en el fenómeno de bronquiolización, donde las cavidades alveolares no ciliadas son invadidas por células basales pulmonares que intentan reparar el daño, sin embargo estás células diferencian a células del epitelio bronquial altamente ciliado, mostrando este curioso enriquecimiento.
Del lado de los gráficos de CC, se realizó para confirmar que estos efectos están enriquecidos principalmente en estructuras que están implicadas contacto y comunicación con la MEC, lo cual es correcto con lo mostrado en los gráficos.

## Complementación con GSEA


Utilizando el objeto del análisis de expresión diferencial con DESeq2 se pudo realizar un ensayo de GSEA. Esto con el fin de obtener un gráfico de dos vías enriquecidas opuestamente en la comparación IPF contra NDC.

**Resultados**: Se graficó a dos vías en un gráfico de GSEA, que se presenta a continuación:

![Gráfico de GSEA. IPF a la izquierda, control a la derecha](https://raw.githubusercontent.com/Hectorjlo/transcriptomica_proyecto_final/refs/heads/main/results/star/GSEA/plots/GSEA_plot.png)

Un gráfico de GSEA nos ofrece la posibilidad de integrar el enriquecimiento de GO para la totalidad de los DEGs sin necesidad de separarlos. Lo que al final da un gráfico mucho más informativo y completo sobre cómo cambia el enriquecimiento de GO desde el IPF hasta el NDC.
En este gráfico, se muestra que la vía GO implicada en la organización de la matriz extracelular (GO:0030198) está altamente enriquecida en el extremo izquierda, que coincide con la condición IPF y a medida que la condición se transforma a NDC este enriquecimiento disminuye, a su vez la vía de la regularización del desarrollo vascular (GO:1901342) tiene un ligero enriquecimiento en la condición IPF, pero en donde sucede realmente su enriquecimiento absoluto más alto es el extremo derecho correspondiente a la condición NDC, esto se debe a que esta vía es necesaria para orquestar la fina estructura vascular alveolar, a medida que la condición se dirige a IPF esta organización se ve afectada lo que se traduce en la pérdida progresiva de la organización, sin embargo el enriquecimiento del lado IPF no desaparece ya que a pesar de su perturbación, el cambio de la estructura de la MEC requiere de la organización de la vascularización.

## Análisis de Interacción proteína-proteína

Para no quedarse solo en un enfoque génico y transcriptómico, se analizó las interacciones de los productos de los transcritos. Para ello se uso STRING (Search Tool For Recurring Instances of Neighbouring Genes) una herramienta para recuperar y mostrar los genes que interactúan con otros repetidamente en el genoma en forma de sus productos proteicos [14]. Se utilizó la totalidad de los DEGs para observar cómo interactúan todos entre sí.

**Resultados:** La herramienta generó una imágen de las interacciones proteína-proteína (IPP) que se muestra a continuación:

![STRING IPP de DEGs. Organizados por clusters](https://raw.githubusercontent.com/Hectorjlo/transcriptomica_proyecto_final/refs/heads/main/results/star/STRING/string_hires_image.png)

A pesar de la alta cantidad de DEGs, los clusters ayudan a identificar grupos definidos que interactúan entre sí. Se observa que la gran cantidad de DEGs interactúa en gran medida entre sí. De STRING se tomó una captura de los primeros doce clusters, y STRING provee en qué mecanismos están implicados:

![STRING detalle de clusterización](https://raw.githubusercontent.com/Hectorjlo/transcriptomica_proyecto_final/refs/heads/main/results/star/STRING/Top_12_clusters_info.jpg)

Nuevamente los resultados vuelvan a apuntar a organización de la MEC, respuesta inmune, quimiocinas resultado de una inflamación crónica, estructuras de tipo ciliadas, deposición de colágeno e incluso eventos de ciclo celular. Este último debido a que en IPF la rigidez, inflamación y desorganización inducen senescencia celular en este tejido.

# Discusión

Actualmente, realizar ensayos puramente biológicos para cuantificar y comparar estadios diferentes tiene la desventaja de que los resultados suelen ser macroscópicos y los posibles blancos terapéuticos quedan ocultos detrás de la complejidad biológica. Por otro lado, ensayos puramente bioinformáticos también carecen de validación en campo y no quedan exentos de los artefactos inherentes a los modelos que se utilizan.
Por ello la integración de estas dos disciplinas es tan importante, combinando las carencias de una son rescatadas por las virtudes de la otra. Debido a esto es fundamentalmente importante comprender la intersección de estas disciplinas, además que la tecnología lo permite y desaprovechar las oportunidades de investigación que ofrecen ambas sería sesgar posibles avances que beneficien a la humanidad.
En este proyecto esta intersección fue probada. La recuperación de los resultados previa fue posible debido a la integración de herramientas bioinformáticas que permiten abstraer la complejidad de simples secuencias de nucleótidos hasta gráficos interpretables de manera biológica. La elección de las herramientas no es trivial, debe ser fundamentada y entendida para que los resultados no sean artefactos propios de la abstracción.
Además por la parte biológica, los resultados obtenidos e interpretados previamente, muestran el gran cambio que implica un proceso de fibrosis en un tejido delicado y altamente organizado como lo es el pulmón. La gran cantidad de DEGs obtenidos apunta justamente como la IPF promueve un cambio drástico en el transcriptoma tisular, mediado por remodeladores de la MEC, respuesta inmune, inflamatoria, desorganización de la vascularización, deposición de fibras de colágeno y reparación desorganizada del tejido, lo que indica que es una enfermedad multifactorial y de terapia complicada dados sus múltiples retroalimentadores.

# Conclusiones 

Se ha de recordar que el uso de las herramientas descritas en el presente documento son esenciales para la obtención de los resultados biológicos interpretados, y su elección debe ser cautelosa para cada caso, los parámetros tienen que tener sustento en el conocimiento biológico para evitar agregar información errónea al final.
En conclusión, la integración de tecnologías bioinformáticas y biológicas es esencial para la elucidación de nuevos conocimientos en las fechas en las que estamos presentes. Es importante notar que no solo es importante la integración a nivel metodológico, sino también a nivel conceptual, entender ambas disciplinas ofrece la conciencia para discernir que lo obtenido no son artefactos bioinformáticos y para descubrir participaciones biológicas antes no consideradas observadas en análisis bioinformáticos.

# Reproducibilidad
Los archivos crudos de lecturas FASTQ fueron descargadas desde archivos SRA provenientes del proyecto [GSE213001](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213001), procesadas en una limpieza por `fastp` v0.24.1, alineadas por `STAR` v2.7.11b contra el genoma de refencia GRCh38.104 obtenido de [GENCODE](https://www.gencodegenes.org/human/release_38.html), los conteos de lecturas fueron cuantificados a genes anotados mediante `featureCounts` v2.0.8, genes con baja expresión fueron eliminados mediante la función de “filterByExpr” proveniente de `edgeR` v4.8.2, para la expresión diferencial se utilizó la librería de `DESeq2` v1.50.2, para el análisis funcional [DAVID](https://davidbioinformatics.nih.gov/summary.jsp) v6.8 fué utilizado, para GSEA se utilizó la librería `clusterProfiler` v4.18.4, y por último para las IPP se usó STRING v12.5.

# Referencias

Babraham Bioinformatics - FastQC A Quality Control tool for High Throughput Sequence Data. (s. f.). https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics, 34(17), i884-i890. https://doi.org/10.1093/bioinformatics/bty560
Dobin, A., Davis, C. A., Schlesinger, F., Drenkow, J., Zaleski, C., Jha, S., Batut, P., Chaisson, M., & Gingeras, T. R. (2012). STAR: ultrafast universal RNA-seq aligner. Bioinformatics, 29(1), 15-21. https://doi.org/10.1093/bioinformatics/bts635
Fu, C., Chen, L., Cheng, Y., Yang, W., Zhu, H., Wu, X., & Cai, B. (2023). Identification of immune biomarkers associated with basement membranes in idiopathic pulmonary fibrosis and their pan-cancer analysis. Frontiers In Genetics, 14, 1114601. https://doi.org/10.3389/fgene.2023.1114601
Huang, D. W., Sherman, B. T., & Lempicki, R. A. (2008). Systematic and integrative analysis of large gene lists using DAVID bioinformatics resources. Nature Protocols, 4(1), 44-57. https://doi.org/10.1038/nprot.2008.211
Huang, G., Gen, Y., Kulur, V., Liu, N., Liu, X., Taghavifar, F., Liang, J., Noble, P. W., & Jiang, D. (2024a). Arrestin beta 1 Regulates Alveolar Progenitor Renewal and Lung Fibrosis. Journal Of Respiratory Biology And Translational Medicine, 1(2), 10006. https://doi.org/10.35534/jrbtm.2024.10006
Huang, G., Gen, Y., Kulur, V., Liu, N., Liu, X., Taghavifar, F., Liang, J., Noble, P. W., & Jiang, D. (2024b). Arrestin beta 1 Regulates Alveolar Progenitor Renewal and Lung Fibrosis. Journal Of Respiratory Biology And Translational Medicine, 1(2), 10006. https://doi.org/10.35534/jrbtm.2024.10006
Huang, X., Yu, W., Tian, J., Zhang, Y., Wei, A., Li, Y., & Chen, S. (2025). Identification and analysis of extracellular matrix and epithelial-mesenchymal transition-related genes in idiopathic pulmonary fibrosis by bioinformatics analysis and experimental validation. Gene, 956, 149464. https://doi.org/10.1016/j.gene.2025.149464
Liao, Y., Smyth, G. K., & Shi, W. (2013). featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics, 30(7), 923-930. https://doi.org/10.1093/bioinformatics/btt656
Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15(12), 550. https://doi.org/10.1186/s13059-014-0550-8
Ma, D., Hua, Y., & Lu, Y. (2024). Complement factor H drives idiopathic pulmonary fibrosis by autocrine C3 regulation, suppressing macrophage phagocytosis and enhancing fibrotic progression. Biochemical And Biophysical Research Communications, 745, 151220. https://doi.org/10.1016/j.bbrc.2024.151220
Ncbi. (s. f.). 08. prefetch and fasterq dump. GitHub. https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump/0681943c0b9e1f2598b2f923bd960cf7626ce235
Pattaroni, C., Begka, C., Cardwell, B., Jaffar, J., Macowan, M., Harris, N. L., Westall, G. P., & Marsland, B. J. (2024). Multi‐omics integration reveals a nonlinear signature that precedes progression of lung fibrosis. Clinical & Translational Immunology, 13(1), e1485. https://doi.org/10.1002/cti2.1485
Sherman, B. T., Hao, M., Qiu, J., Jiao, X., Baseler, M. W., Lane, H. C., Imamichi, T., & Chang, W. (2022). DAVID: a web server for functional enrichment analysis and functional annotation of gene lists (2021 update). Nucleic Acids Research, 50(W1), W216-W221. https://doi.org/10.1093/nar/gkac194
Snel, B. (2000). STRING: a web-server to retrieve and display the repeatedly occurring neighbourhood of a gene. Nucleic Acids Research, 28(18), 3442-3444. https://doi.org/10.1093/nar/28.18.3442


# Anexo 
n.1:
```bash
cat <(awk -v OFS='\t' -v FPAT=' *"[^"]*" *|[^,]*' 'NR == 1 {print $1, $22, $33, $15, $3, $25, $48}' SraRunTable.csv) <(awk -v OFS='\t' -v FPAT=' *"[^"]*" *|[^,]*' '($15 == "IPF" || $15 == "NDC") && (($15 == "IPF" && $33 == "Base") || $15 == "NDC") && !seen[$22]++ {print $1, $22, $33, $15, $3, $25, $48}' SraRunTable.csv | sort -k4) > extracted_samples.tsv

awk 'NR > 1 {print $1}' extracted_samples.tsv > SRR_list.txt
```
n.2:
```bash

join <(for file in *.fastq; do
echo "${file%%.fastq}"
	done | sort) <(awk -v OFS="" 'NR > 1 {print $1, "\t", $4}' ../../data/GSE213001/extracted_samples.tsv | sort -k1) | awk -v OFS="" '{print "mv ", $1, ".fastq", " ", $2, "_", $1, ".fastq"}' | bash
``` 
n.3:
```bash
STAR --runMode genomeGenerate \
     --genomeDir STAR_index \
     --genomeFastaFiles GRCh38.p13.genome.fa \
     --sjdbGTFfile gencode.v38.chr_patch_hapl_scaff.annotation.gtf \
     --sjdbOverhang 87 \
     --runThreadN 40

```
n.4:
```bash
src/auto_BAMic.sh results/star/alignments results/star/alignments &> results/star/alignments/BAMs.log

rm -rf results/star/alignments/*.sam
```

n.5:
```bash
awk -v OFS="\t" 'NR > 1  {print}' results/star/feature_counts/ft_counts.tsv | cut -f1,7-40 > results/star/feature_counts/counts_matrix.tsv

```

n.6:
```bash
cat data/GENCODE_GRCh38.p13_104/gencode.v38.chr_patch_hapl_scaff.annotation.gtf | awk '$3 == "gene"' | grep -oP 'gene_id "\K[^"]+|gene_name "\K[^"]+' | paste -d'\t' - - > data/GENCODE_GRCh38.p13_104/gene_id.gene_name.txt
```
