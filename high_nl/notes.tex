
CALCULO DE ENERGÍA CON SHIFTS:
Hice cálculos con dos shifts diferentes:
1) shiftls: el shift se hace sobre los términos (en carpeta ishftls/).
2) shiftic: el shift se hace sobre los niveles (en carpeta ishftic/).
Los resultados los obtuve usando el flag ISHFTLS/IC>1, lo que significa
que el programa adopta las energías observadas (dadas en el archivo 
SHFTLS/IC) relativas al nivel fundamental, y luego el programa prosigue 
a hacer iteraciones con las correcciones (ISHFTLS o ISHFTIC) sobre H(IC). 
Las iteraciones ISHTFLS se implementan como correcciones a las energías
de los términos (TEC) al hamiltoniano IC.

Los mejores resultados (Aki & Ei) los obtuve haciendo un shift en LS
(ver comparación en "compare_shifts.ipynb"). La desviación encontrada
con ISHFTIC>1 fue de ~0.2\% (tanto para los niveles como para los términos),
mientras que con ISHFTLS>1, la desviación ~1E-04\%.

ESCRITURA DE SHFTLS/IC:
La escritura de los archivos SHFTLS/IC es automático y se hace con 
"write_ishft.ipynb". Este programa tiene como inputs los archivos
"NIST_cfgs.dat" y "NIST_levels.dat". El formato de estos archivos es
el siguiente: 
>> NIST_cfgs.dat: 
  "i CFG"
  > i   - índice que utiliza AS para cada configuración (según input das). 
  > CFG - configuracion electronica con formato NIST
>> NIST_levels.dat: 
  "Configuration Term J Level(Ry)"
  > Configuration - configuración electronica con formato NIST
  > Term          - término espectroscópico
  > J             - duh
  > Level(Ry)     - energía del nivel en unidades de Rydberg
Ambos archivos tiene un header de 3 líneas.

SHIFT DE LEVELS/TERMS SIN DATOS DE NIST:
Un problema con el que me encontré al hacer los shifts fue la falta de 
algunos valores experimentales de energía. Entonces, cuando hacía los 
shifts, lo que occuría es que los términos/niveles que no tenían valores 
observados tomaban los valores del cálculo *sin* shift y quedaban muy 
por debajo. Finalmente, lo que hice para corregir este problema fue
hacer un promedio de las desviaciones (entre los valores de AS y NIST)
de los niveles más altos y crear shifts artificiales para los niveles 
sin observaciones (ver "write_ishft.ipynb"). Por supuesto, los shifts 
artificiales quedaron montados sobre los niveles teóricos.

CALCULO DE TRANSICIONES RADIATIVAS USANDO SHIFTLS
Hice dos cálculo de transiciones radiativas. 
1) Consideré sólo las transiciones eléctricas dipolares (E1). Hice la 
lectura de dos archivos de escritura de autostructure y encontré una 
diferencia en el número de transiciones impresos. En "olg" AS escribe 
todas las transciones, mientras que en "oic" AS imprime transicones tal 
que Aki>0.01. Finalmente opté por tomar sólo las transiciones del archivo
"oic".
2) Consideré las transiciones E1, E2/M1 y E3/M2. Nuevamente encontré una
diferencia en el número de transciones impresas en los archivos "olg" y
"oic". 

