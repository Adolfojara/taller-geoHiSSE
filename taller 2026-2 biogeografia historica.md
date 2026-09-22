# Tutorial: Uso de GeoHiSSE con árboles obtenidos desde TimeTree

**Autor:** (tu nombre)  
**Fecha:** (fecha)  
**Licencia:** CC BY 4.0

---

## Índice

1. [Instalación y carga de paquetes](#paso-1-instalación-y-carga-de-paquetes)
2. [Obtener el árbol desde TimeTree](#paso-2-obtener-el-árbol-desde-timetree)
3. [Cargar y verificar el árbol en R](#paso-3-cargar-y-verificar-el-árbol-en-r)
4. [Preparar los datos de rangos geográficos](#paso-4-preparar-los-datos-de-rangos-geográficos)
5. [Construir la matriz de transición](#paso-5-construir-la-matriz-de-transición)
6. [Especificar tasas de diversificación](#paso-6-especificar-tasas-de-diversificación)
7. [Ejecutar GeoHiSSE](#paso-7-ejecutar-geohisse)
8. [Comparar modelos](#paso-8-comparar-modelos)
9. [Interpretar resultados](#paso-9-interpretar-resultados)
10. [Visualización y análisis adicionales](#paso-10-visualización-y-análisis-adicionales)
11. [Resumen del flujo](#resumen-del-flujo)
12. [Limitaciones y advertencias sobre TimeTree](#limitaciones-y-advertencias-sobre-timetree)

---

## Introducción

**GeoHiSSE** (Geographic Hidden State Speciation and Extinction) es un modelo de Markov oculto implementado en el paquete `hisse` de R. Permite investigar cómo la distribución geográfica de los linajes influye en sus tasas de diversificación, incorporando **estados ocultos** que capturan heterogeneidad no observada en las tasas, tanto entre regiones como dentro de ellas.

En este tutorial, el árbol filogenético será obtenido por los estudiantes desde **TimeTree** (www.timetree.org), una base de datos pública que integra miles de estudios de datación molecular y ofrece árboles calibrados en el tiempo en formato **Newick**.

**Requisitos previos:**

- R ≥ 4.0 y RStudio
- Conexión a internet (para instalar paquetes y descargar el árbol)
- Conocimientos básicos de R
- Datos de rangos geográficos codificados para las especies de interés

---

## Paso 1: Instalación y carga de paquetes

**Qué se hace:** Instalar la versión de desarrollo de `hisse` desde GitHub y cargar `hisse`, `diversitree`, `ape` y `phytools`.

**Por qué:** `hisse` contiene las funciones centrales de GeoHiSSE; `diversitree` aporta infraestructura para verosimilitud; `ape` y `phytools` se usan para leer y manipular árboles Newick.

```r
# Instalar devtools si no lo tienes
install.packages("devtools")

# Instalar hisse desde GitHub
library(devtools)
install_github(repo = "thej022214/hisse", ref = "master")

# Cargar los paquetes necesarios
suppressWarnings(library(hisse))
suppressWarnings(library(diversitree))
library(ape)
library(phytools)
```

Al cargar `hisse` se cargan automáticamente dependencias como `ape`, `deSolve`, `GenSA`, `subplex` y `nloptr`.

---

## Paso 2: Obtener el árbol desde TimeTree

**Qué se hace:** Buscar el grupo de interés en TimeTree y descargar el árbol en formato Newick.

**Por qué:** TimeTree ofrece árboles ya calibrados en millones de años (Ma), lo que evita tener que hacer una ultrametricación "a ciegas". Sin embargo, hay que elegir bien el tipo de árbol y verificar su calidad antes de usarlo.

### 2.1 Pasos en el sitio web

1. Visitar [www.timetree.org](http://www.timetree.org)
2. Buscar el grupo objetivo (especie, género o clado mayor). TimeTree acepta nombres científicos y comunes.
3. Seleccionar el estudio que contenga todos los taxones de interés.
4. Descargar el árbol en formato Newick (`.nwk`).

### 2.2 Elegir el tipo correcto de árbol

TimeTree ofrece tres versiones para cada consulta:

| Tipo | Qué hace | ¿Adecuado para GeoHiSSE? |
|------|----------|---------------------------|
| **Unsmoothed (sin suavizar)** | Conserva las politomías originales sin resolverlas | **Sí, el más recomendado** |
| **Smoothed (suavizado)** | Resuelve aleatoriamente las politomías en dicotomías | Aceptable, pero introduce topología arbitraria |
| **Interpolated (interpolado)** | Inserta especies sin datos según su posición taxonómica | **No recomendado**: distorsiona longitudes de rama |

**Recomendación:** descargar siempre la versión **unsmoothed**.

---

## Paso 3: Cargar y verificar el árbol en R

**Qué se hace:** Leer el archivo Newick en R y comprobar sus propiedades básicas.

**Por qué:** Antes de continuar, hay que asegurarse de que el árbol está enraizado, es ultramétrico (o casi) y no contiene errores de formato.

```r
# Cargar el árbol de TimeTree
phy <- read.tree("timetree_mi_grupo.nwk")

# Verificar propiedades
print(phy)                 # resumen general
phy$tip.label              # nombres de las especies
is.rooted(phy)             # debe ser TRUE
is.binary(phy)             # FALSE es normal en árboles unsmoothed
is.ultrametric(phy)        # se espera TRUE
```

### 3.1 Si `is.ultrametric()` devuelve FALSE

A veces el árbol es ultramétrico en teoría, pero los errores de precisión numérica hacen que las distancias raíz-hoja difieran en valores muy pequeños.

```r
# Verificar magnitud del problema
distancias <- dist.nodes(phy)[1:length(phy$tip.label), 1]
range(distancias)   # si el rango es < 1e-6, es solo precisión numérica

# Corregir extendiendo las ramas terminales
phy <- force.ultrametric(phy, method = "extend")
is.ultrametric(phy)  # ahora TRUE
```

**Importante:** `force.ultrametric(method = "extend")` solo debe usarse como corrección numérica. **No** es un método de calibración temporal.

### 3.2 Sobre las politomías

Los árboles unsmoothed de TimeTree suelen contener politomías (nodos con más de dos descendientes). GeoHiSSE las maneja bien; **no es necesario** convertir el árbol a binario con `multi2di()`, ya que eso introduciría topología arbitraria.

---

## Paso 4: Preparar los datos de rangos geográficos

**Qué se hace:** Construir un data frame con dos columnas: nombre de especie y código de rango geográfico.

**Por qué:** GeoHiSSE requiere que cada especie esté codificada en un modelo de dos áreas: `0` = rango amplio (ambas áreas), `1` = solo área 0, `2` = solo área 1.

```r
# Estructura esperada
sim.dat <- data.frame(
  taxon  = c("Panthera_leo", "Panthera_tigris", "Panthera_pardus"),
  ranges = c(1, 1, 0)   # 0 = amplio, 1 = área 0, 2 = área 1
)
```

### 4.1 Sincronizar nombres con el árbol

TimeTree usa típicamente el formato `Género_especie` con guion bajo. Si tu tabla usa espacios o abreviaturas, hay que homogeneizar.

```r
# Normalizar nombres
sim.dat$taxon <- gsub(" ", "_", sim.dat$taxon)

# Verificar coincidencias
matched <- sim.dat$taxon %in% phy$tip.label
cat("Especies coincidentes:", sum(matched), "de", nrow(sim.dat), "\n")

# Ver qué especies no coinciden
sim.dat$taxon[!matched]
```

Si más del 20% de las especies no coinciden, revisa la taxonomía antes de continuar. GeoHiSSE ignorará silenciosamente las especies no presentes en el árbol.

### 4.2 Reordenar los datos para que coincidan con el árbol

```r
# Opcional: ordenar sim.dat según phy$tip.label
sim.dat <- sim.dat[match(phy$tip.label, sim.dat$taxon), ]
```

### 4.3 Calcular la fracción de muestreo `f`

Si el árbol **no contiene todas las especies vivas** del grupo, hay que estimar la fracción de muestreo para cada estado geográfico.

```r
# Ejemplo concreto:
# Especies de rango amplio: 12 en árbol / 20 totales → 0.60
# Especies solo del área 0: 25 en árbol / 40 totales → 0.625
# Especies solo del área 1: 23 en árbol / 40 totales → 0.575
f <- c(0.60, 0.625, 0.575)
```

Si **todas** las especies vivas están en el árbol, entonces `f = c(1, 1, 1)`.

---

## Paso 5: Construir la matriz de transición

**Qué se hace:** Usar `TransMatMakerGeoHiSSE()` para generar los índices de la matriz de transición geográfica y oculta.

**Por qué:** Esta matriz define qué transiciones están permitidas y qué parámetros se comparten entre estados.

```r
# Modelo GeoSSE sin estados ocultos
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits = 0)

# Modelo GeoHiSSE con 1 estado oculto
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits = 1)

# Modelo nulo (transiciones simplificadas)
trans.rate.null <- TransMatMakerGeoHiSSE(hidden.traits = 1, make.null = TRUE)
```

**Argumentos relevantes:**

- `hidden.traits`: número de estados ocultos (0 para GeoSSE clásico).
- `make.null`: iguala tasas entre estados ocultos.
- `include.jumps`: permite transiciones directas 0↔1 sin pasar por 01.
- `separate.extirpation`: distingue la extirpación desde 01 hacia 0 vs. hacia 1.

La matriz devuelta es de **índices**: números iguales comparten un mismo parámetro, `NA` indica transición prohibida, `0` indica parámetro fijado.

---

## Paso 6: Especificar tasas de diversificación

**Qué se hace:** Definir los vectores `turnover` y `eps` (extinction fraction).

**Por qué:** Los números en estos vectores son **índices**, no valores. Índices iguales implican que esos parámetros se estiman como uno solo.

```r
# Sin estados ocultos
turnover <- c(1, 2, 3)   # s0, s1, s01
eps      <- c(1, 2)      # x0, x1

# Con 1 estado oculto, todos libres
turnover <- c(1, 2, 3, 4, 5, 6)   # s0A, s1A, s01A, s0B, s1B, s01B
eps      <- c(1, 2, 3, 4)         # x0A, x1A, x0B, x1B

# Modelo restringido (tasas iguales entre estados ocultos)
turnover <- c(1, 2, 3, 1, 2, 3)
eps      <- c(1, 2, 1, 2)
```

**Regla:** la longitud del vector debe coincidir con el número de parámetros del modelo. Con 1 estado oculto: `turnover` de longitud 6, `eps` de longitud 4.

---

## Paso 7: Ejecutar GeoHiSSE

**Qué se hace:** Llamar a `GeoHiSSE()` para ajustar el modelo por máxima verosimilitud.

**Por qué:** Es el paso central de estimación.

```r
modelo <- GeoHiSSE(
  phy            = phy,
  data           = sim.dat,
  f              = f,                    # fracciones calculadas en el Paso 4.3
  turnover       = c(1, 2, 3, 4, 5, 6),
  eps            = c(1, 2, 3, 4),
  hidden.states  = TRUE,
  trans.rate     = trans.rate,
  sann           = TRUE,
  sann.its       = 1000
)
```

**Parámetros clave:**

- `f`: fracciones de muestreo por estado geográfico.
- `hidden.states`: incluir o no estados ocultos.
- `sann`: usar *simulated annealing* para búsqueda global (recomendado).
- `assume.cladogenetic`: si se asume especiación por vicarianza (por defecto `TRUE`).
- `root.type`: manejo del estado en la raíz (`"madfitz"` o `"herr_als"`).

**Tiempo de cómputo:** GeoHiSSE puede tardar horas o días. Para probar la configuración, usa primero `sann.its = 100`.

---

## Paso 8: Comparar modelos

**Qué se hace:** Ajustar varios modelos con hipótesis distintas y compararlos con AIC/AICc.

**Por qué:** Un AIC aislado no dice nada; hay que comparar hipótesis competidoras.

```r
# Modelo 1: GeoSSE sin estados ocultos
mod1 <- GeoHiSSE(phy, sim.dat, f = f,
                 turnover = c(1,2,3), eps = c(1,2),
                 hidden.states = FALSE, trans.rate = trans.rate.null)

# Modelo 2: GeoHiSSE con 1 estado oculto
mod2 <- GeoHiSSE(phy, sim.dat, f = f,
                 turnover = c(1,2,3,4,5,6), eps = c(1,2,3,4),
                 hidden.states = TRUE, trans.rate = trans.rate)

# Modelo 3: independiente de región (s01 removido)
mod3 <- GeoHiSSE(phy, sim.dat, f = f,
                 turnover = c(1,1,0,2,2,0),   # s01 fijado a 0
                 eps      = c(1,1,2,2),
                 hidden.states = TRUE, trans.rate = trans.rate)

# Comparación
AIC(mod1, mod2, mod3)
```

**Principios:**

- En el modelo independiente de región, `s01` debe fijarse a `0`.
- Usar **pesos AIC** además de las diferencias.
- Evitar modelos con demasiados parámetros libres (>~15) para reducir falsos positivos.

---

## Paso 9: Interpretar resultados

**Qué se hace:** Extraer tasas y convertirlas a especiación/extinción.

**Por qué:** GeoHiSSE estima turnover (τ = λ + μ) y extinction fraction (ε = μ/λ), no λ y μ directamente.

```r
summary(modelo)
modelo$solution
modelo$AIC
modelo$AICc
```

**Conversión:**

- λ = τ / (1 + ε)
- μ = τ × ε / (1 + ε)

**Unidades:** como el árbol de TimeTree está calibrado en millones de años, las tasas tienen unidades de **Ma⁻¹**. Por ejemplo, λ = 0.1 Ma⁻¹ significa 0.1 eventos de especiación por linaje por millón de años.

**Interpretación de estados ocultos:** representan fuentes de variación en las tasas **no relacionadas con el rango geográfico**. Si el modelo los favorece, hay heterogeneidad dentro de cada región.

---

## Paso 10: Visualización y análisis adicionales

```r
# Reconstrucción de rangos ancestrales
recon <- MarginReconGeoSSE(modelo, phy = phy, data = sim.dat)

# Graficar
plot(recon)
```

También existe `GetModelAveRates()` para promediar tasas entre modelos cuando ninguno domina claramente.

---

## Resumen del flujo

| Paso | Acción | Propósito |
|------|--------|-----------|
| 1 | Instalar y cargar paquetes | Preparar entorno |
| 2 | Descargar árbol de TimeTree (unsmoothed) | Obtener filogenia calibrada |
| 3 | Verificar y corregir el árbol en R | Asegurar ultrametricidad |
| 4 | Preparar datos de rango y calcular `f` | Alinear datos con el árbol |
| 5 | Construir matriz de transición | Definir transiciones permitidas |
| 6 | Especificar `turnover` y `eps` | Definir restricciones |
| 7 | Ejecutar `GeoHiSSE()` | Estimar parámetros |
| 8 | Comparar modelos con AIC | Seleccionar mejor hipótesis |
| 9 | Interpretar tasas | Traducir a biología |
| 10 | Visualizar y reconstruir | Analizar resultados |

---

## Limitaciones y advertencias sobre TimeTree

1. **Cobertura incompleta:** muchos grupos carecen de datos moleculares suficientes; TimeTree puede devolver árboles con pocos representantes.
2. **Los valores de soporte no son bootstrap:** TimeTree reporta la fracción de estudios que apoyan un nodo, no la confianza filogenética clásica. GeoHiSSE no usa estos valores directamente, pero la fiabilidad del árbol sí importa.
3. **Politomías frecuentes:** en árboles unsmoothed, esto es normal y manejable.
4. **Fracción de muestreo:** casi siempre menor que 1, lo que obliga a estimar `f` con cuidado.
5. **Si el árbol es smoothed:** la topología es arbitraria en los nodos politómicos originales; hay que reportar esto como una fuente de incertidumbre.

---

## Conclusión

Usar TimeTree como fuente del árbol tiene ventajas claras: es rápido, reproducible y ofrece árboles calibrados en tiempo absoluto, lo que permite interpretar las tasas de GeoHiSSE en Ma⁻¹. Sin embargo, exige tres verificaciones obligatorias: (1) elegir la versión *unsmoothed*, (2) comprobar la ultrametricidad y corregirla si es solo un problema numérico, y (3) calcular correctamente la fracción de muestreo `f`. Con estas precauciones, el flujo completo —desde la instalación de `hisse` hasta la comparación de modelos y la interpretación biológica— puede ejecutarse íntegramente en R y producir resultados defendibles.

---

## Referencias

- Beaulieu, J. M., & O'Meara, B. C. (2016). Detecting hidden diversification shifts in models of trait-dependent speciation and extinction. *Systematic Biology*, 65(4), 583–601.
- Kumar, S., Stecher, G., Suleski, M., & Hedges, S. B. (2017). TimeTree: a resource for timelines, timetrees, and divergence times. *Molecular Biology and Evolution*, 34(7), 1812–1819.
- Nakov, T., Beaulieu, J. M., & Alverson, A. J. (2019). Insights into global planktonic diatom diversity: an integrated analysis of fossil and molecular data. *Frontiers in Marine Science*, 6, 329.
- Documentación del paquete `hisse`: https://github.com/thej022214/hisse
