# TFM_IBD_GSE214695

Pipeline de análisis de single-cell RNA-seq (scRNA-seq) en mucosa colónica de enfermedad inflamatoria intestinal (EII), a partir del dataset público **GSE214695**, comparando controles sanos (HC), colitis ulcerosa (UC) y enfermedad de Crohn (CD). Además, expresión diferencial pseudobulk y modelo de clasificación por machine learning.

Trabajo Fin de Máster — Paula Damián.

## Documento del TFM

La justificación metodológica, la discusión biológica de los resultados y las conclusiones se encuentran en la memoria del TFM (documento Word), no en este repositorio. Este README documenta el aspecto técnico.

## Origen de los datos

GSE214695 en GEO, 18 muestras (6 HC, 6 UC, 6 CD), matrices de conteos crudas por muestra más un archivo de anotación celular oficial del propio estudio.
| Publicación original | Garrido-Trigo et al., *Nat Commun* (2023), PMID [37495570](https://pubmed.ncbi.nlm.nih.gov/37495570/) |

Los datos no se distribuyen en este repositorio. Se descargan automáticamente al ejecutar `notebooks/00_download_data.ipynb`.

## Estructura del repositorio

Todo el código y la configuración viven en este repositorio (TFM_IBD_GSE214695/). Los datos crudos y procesados están en una carpeta aparte en Google Drive (IBD_TFM/), que no está versionada. Cada notebook resuelve las rutas combinando dos variables:

- REPO_ROOT: dónde está clonado este repo. Se usa para leer config/params.yaml y para guardar figuras y tablas (reports/).
- PROJECT_ROOT: dónde están los datos (IBD_TFM). Se usa para todo lo que sea un objeto .h5ad, crudo o intermedio.

```
TFM_IBD_GSE214695/
├── config/
│   ├── params.yaml            
│   └── qc_thresholds.yaml     
├── data/                        
│   ├── raw/
│   └── interim/
│       ├── 00_raw_h5ad/
│       ├── 01_qc/filtered/
│       ├── 02_normalized/
│       ├── 03_diagnostics/
│       ├── 04_integrated/
│       ├── 05_clustered/
│       ├── 06_cluster_diagnostics/
│       ├── 07_annotated/
        ├──  08_deg_pseudobulk/
│       └── 09_ml_classification/
├── notebooks/
│   ├── 00_download_data.ipynb
│   ├── 01_qc_interactive.ipynb
│   ├── 02_normalization.ipynb
│   ├── 03_diagnostics.ipynb
│   ├── 04_integration.ipynb
│   ├── 05_clustering.ipynb
│   ├── 06_cluster_diagnostics.ipynb
│   ├── 07_annotation.ipynb
│   ├── 08_deg_pseudobulk.ipynb
│   └── 09_ml_classification.ipynb
├── scripts/
│   └── 01_qc.py                
├── reports/
│   ├── figures/                
│   └── tables/
├── checks/
│   └── 0X_anexo_hvg.ipynb
├── environment.yml
├── environment_deg_ml.yaml
├── CITATION.cff
├── LICENSE
└── README.md
```

## Flujo del pipeline

```
00_download_data → 01_qc → 02_normalization → 03_diagnostics
    → 04_integration → 05_clustering → 06_cluster_diagnostics → 07_annotation
    → 08_deg_pseudobulk → 09_ml_classification
```

Paso	Qué hace
00_download_data	Descarga GSE214695 de GEO y convierte cada muestra a .h5ad
01_qc	Control de calidad y filtrado, umbrales fijados por muestra
02_normalization	Concatena las 18 muestras, normaliza, HVGs, PCA
03_diagnostics	Añade la condición (HC/UC/CD) y diagnostica visualmente el batch effect antes de integrar
04_integration	Corrige el batch effect con Harmony
05_clustering	Leiden a varias resoluciones, se fija la de 0.8
06_cluster_diagnostics	Revisión de los clusters: dominancia por muestra, QC conjunto, identidad por marcadores, tabla de decisión (qué cluster se queda, cuál se fusiona o se descarta)
07_annotation	Anotación celular final, cruzada con la anotación oficial del estudio
08_deg_pseudobulk Agregación pseudobulk por paciente y tipo celular, expresión diferencial (DESeq2), tabla de features para ML 09_ml_classification Clasificación HC/UC/CD (LOOCV anidado), test de permutación, interpretabilidad (SHAP)

01_qc existe en dos versiones que hacen lo mismo por caminos distintos:`config/qc_thresholds.yaml` contiene los umbrales de QC ya decididos y congelados para las 18 muestras. `scripts/01_qc.py` los aplica de forma determinista y reproducible. `notebooks/01_qc_interactive.ipynb` permite explorar otros umbrales de forma interactiva sin afectar al resultado oficial del pipeline.

## Entornos de ejecución

Este proyecto usa dos entornos Conda independientes, uno por fase del pipeline:

Notebooks	Entorno	Contenido
00 a 07 (descarga → anotación celular)	environment.yml	Scanpy, AnnData, Harmony, Leiden — ecosistema de análisis de células individuales
08 a 09 (DEG pseudobulk → ML)	environment_deg_ml.yaml	PyDESeq2, formulaic, scikit-learn — modelado estadístico y clasificación

Cada notebook instala su propio entorno en la primera celda (vía condacolab + conda env update), así que no hace falta activar nada manualmente fuera de Colab: basta con ejecutar los notebooks en orden dentro del bloque correspondiente.

## Instalación

**Entorno local**

```bash
conda env create -f environment.yml
conda activate tfm-ibd-scrnaseq
```

**Google Colab**

Colab trae preinstalados algunos paquetes (pytensor, google-colab) que entran en conflicto con las versiones que necesita este pipeline. Por eso cada notebook empieza con una celda que monta el entorno conda y reinicia el kernel (se mostrará "Your session crashed for an unknown reason" justo después, no es un fallo).

Para lanzar `scripts/01_qc.py` desde Colab hace falta lo mismo:
```
# 1)  montar Drive
from google.colab import drive
drive.mount('/content/drive')

!pip install -q condacolab
import condacolab
condacolab.install()
```
```
# 2) Instalación del entorno
!conda env update -n base -f /content/drive/MyDrive/TFM_IBD_GSE214695/environment.yml -q
```
```
!python "/content/drive/MyDrive/TFM_IBD_GSE214695/scripts/01_qc.py"
```

Todas las versiones exactas están fijadas en `environment.yml` y se repiten, idénticas, en la celda de instalación de cada notebook.

## Ejecución del proyecto

1. Clonar el repositorio.
2. Tener la carpeta de datos accesible en la ruta que indica `TFM_ROOT` (variable de entorno) o el valor por defecto de `config/params.yaml`.
3. Ejecutar los notebooks en orden.
4. Para reproducir el QC de forma automática, ejecutar `scripts/01_qc.py` en lugar del notebook interactivo.

## Reproducibilidad

- Rutas y parámetros centralizados en `config/params.yaml`.
- Umbrales de QC congelados en `config/qc_thresholds.yaml`.
- Versiones de dependencias fijadas en `environment.yml`.
- Datos de origen públicos (GEO), descargados automáticamente.

## Licencia y cita

Distribuido bajo licencia [MIT](LICENSE). Para citar este trabajo, ver [CITATION.cff](CITATION.cff).