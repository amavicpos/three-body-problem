> Resolución del problema de los tres cuerpos a partir del sistema
> Kepler-16

En este trabajo se ha resuelto el problema de los tres cuerpos, mediante
el uso de la herramienta Matlab, con los métodos de Euler y Verlet. Se
ha resuelto la órbita del Sistema Kepler-16 con

> exactitud, con un error relativo 𝐸𝑟𝑟𝑜𝑟% = (4.92 × 10−8)% en la energía
> mecánica.
>
> Se han buscado nuevas órbitas adaptando las condiciones de inicio. Con
> esto, se han obtenido resultados de distinta complejidad y analizado
> los datos obtenidos sobre ellos. Además, se ha observado que la
> presencia de órbitas no estables es mucho mayor, debido a que leves
> cambios en las condiciones de inicio afectan a la órbita en gran
> medida. Dicho esto, se ha comprobado que los sistemas binarios
> albergando planetas, que recorren órbitas tanto internas como
> externas, son teóricamente posibles.
>
> En resumen, este acercamiento numérico ofrece resultados del problema
> con exactitud, además de mostrar aspectos como la conservación de la
> energía mecánica y momento angular.
>
> Autores: Martín San Juan y Amaia Vicario
>
> Filiación: Universidad Autónoma de Madrid
>
> Fecha: 4 de mayo del 2021

1

> 1.[Introducción:]{.underline}

El problema de los tres cuerpos consiste en resolver las posiciones de
tres cuerpos libres, en movimiento debido a la atracción gravitatoria
que ejerce cada uno. Para ello, se ha usado la ecuación de Newton:

+-----------------+-----------------+-----------------+-----------------+
| 𝑑2𝑟𝑖            | 𝑚𝑗𝑟𝑖𝑗           |                 | > 𝑚𝑘𝑟𝑖𝑘 (1)     |
+=================+=================+=================+=================+
| 𝑑𝑡2 = 𝐺(        | 3\              | \+              | > 3)\           |
|                 | 𝑑𝑖𝑗             |                 | > 𝑑𝑖𝑘           |
+-----------------+-----------------+-----------------+-----------------+

Donde 𝑑𝑡2 es la aceleración del objeto, 𝑟⃑ su posición, 𝑚𝑗y 𝑚𝑘 son las
masas de los otros

cuerpos, 𝐺es la constante de gravitación universal; 𝑑𝑖𝑗 y 𝑑𝑖𝑘 son las
distancias entre el cuerpo y los otros dos; y 𝑟⃑⃑⃑ y 𝑟𝑖𝑘⃑⃑⃑⃑ son los vectores
que unen al objeto con los otros dos cuerpos.

El interés que tiene este problema es que no tiene una solución
analítica general. Exceptuando algunas soluciones encontradas, la
mayoría de los casos presentan órbitas no periódicas o inestables.

En este trabajo se parte del cálculo de la órbita del sistema Kepler-16;
después, se varían las condiciones iniciales de las que parten los tres
cuerpos para observar posibles soluciones de este problema y estudiar su
estabilidad.

> 2.[Ecuaciones de movimiento]{.underline}:

Como se ha comentado, se necesita la ecuación de Newton (1) para obtener
las aceleraciones de cada cuerpo, y así resolver sus posiciones.

Aparte de eso, dado que los cuerpos se aceleran, el sistema se desplaza
a medida que avanza el tiempo. Por esta razón, se han calculado las
posiciones y energías con las distancias respecto al centro de masas,
que constituye el sistema de referencia en el que el momento lineal
total es cero. Las fórmulas para obtener las posiciones relativas al
centro de masas son las siguientes:

+-----------------------+-----------------------+-----------------------+
|                       | 3\                    | \(2\)                 |
|                       | ∑𝑚𝑖∙𝑟𝑖                |                       |
+=======================+=======================+=======================+
| 𝑟𝐶𝑀=                  | ∑𝑚𝑖3                  | > ; 𝑟𝑟𝑒𝑙= 𝑟𝑖−𝑟𝐶𝑀      |
|                       |                       |                       |
|                       | 𝑖                     |                       |
+-----------------------+-----------------------+-----------------------+

Donde 𝑟𝐶𝑀 es la posición del centro de masas y 𝑟𝑟𝑒𝑙 su posición
relativa. Además, se ha estudiado la estabilidad de las órbitas
observando la conservación de la energía total y del momento angular. La
energía total (E) es la suma de las energías cinética (*K*) y potencial
(*U*) totales, descritas por las fórmulas:

+-----------------------+-----------------------+-----------------------+
| > 3\                  | > 3 (3) ; 𝑈𝑖=         | \(4\)                 |
| > 𝐾𝑖= 1 2 𝑚𝑖∙𝑣(𝑟)⃑⃑⃑⃑⃑⃑⃑⃑ 𝑖2  | > 𝐺𝑚𝑖𝑚𝑗\|𝑟𝑖−𝑟𝑗\| +    |                       |
| > ; 𝐾𝑡𝑜𝑡= ∑𝐾𝑖         | > 𝐺𝑚𝑖𝑚𝑘𝑘\| ; 𝑈𝑡𝑜𝑡=    |                       |
|                       | > ∑𝑈𝑖𝑖                |                       |
| 𝑖                     |                       |                       |
+=======================+=======================+=======================+
+-----------------------+-----------------------+-----------------------+

El momento angular es en esencia un producto vectorial que describimos
como:

> 𝐿𝑖= 𝑟𝑖× (𝑚𝑖∙𝑣𝑖); 𝐿𝑡𝑜𝑡= \|𝐿⃑𝑖+ 𝐿⃑𝑗+𝐿⃑𝑘\|(5)

2

Siendo 𝐿𝑡𝑜𝑡la norma de la suma de los momentos de cada cuerpo. La
dirección de *Li* viene determinada por la regla de la mano derecha.

> 3.[Métodos y desarrollo del código]{.underline}:

En este trabajo se ha hecho uso del programa Matlab para calcular las
órbitas numéricamente; es decir, creando un algoritmo que simule una
aproximación de la solución buscada. Con el fin de calcular las
posiciones y velocidades de los cuerpos en cada instante, se comenzó
aplicando tanto el método Euler\[5\] como el método Verlet\[6\] de
velocidades.

Como se observa de los cálculos hechos con distintos métodos, la
diferencia de tiempo estimada en el cálculo es menor para Euler:
𝑡𝑉𝑒𝑟𝑙𝑒𝑡= 15.46 𝑠 y 𝑡𝑒𝑢𝑙𝑒𝑟= 11.29 𝑠. Sin embargo, el método Verlet
muestra una mejor aproximación a la solución1; Verlet muestra un error
relativo de 𝑉𝑒𝑟𝑙𝑒𝑡 𝐸𝑟𝑟𝑜𝑟 % = (4.92 × 10−8)% en la energía mecánica,
mientras que 𝐸𝑢𝑙𝑒𝑟 𝐸𝑟𝑟𝑜𝑟 % = (6.32 × 10−4)%, valor considerablemente
mayor. Por esta razón, se usa Verlet para las órbitas finales,
exceptuando el caso de las figuras con el segundo centro de masas, las
cuales se calculan con Euler.

La razón para dejarlas así es que el código inicial se escribió con
Euler, y los resultados eran prácticamente iguales que con Verlet, por
lo que, por razones de comodidad (cambiar el código para volver a añadir
los centros de masas era bastante complicado), se han dejado con Euler.

Con el fin de optimizar el cálculo y que sea lo más sencillo posible, se
agrupan las posiciones y velocidades en matrices de dimensión 3xNS (NS:
número de pasos), una para cada componente del parámetro (*x,y,z*),
donde cada columna indica la posición o velocidad de los tres planetas
en cada paso y cada fila indica la evolución de un planeta respecto a un
parámetro. Para las masas de los planetas se utiliza un vector fila, de
manera que las operaciones matriciales cuadren. Agrupar las masas
minimiza el número de veces que se ejecuta Verlet y por lo tanto hace el
código más eficiente.

Para obtener las órbitas de los cuerpos son necesarias ciertas
condiciones iniciales, su periodo orbital y su radio, así como su masa.
La velocidad se obtiene por aproximación a órbitas

circulares (𝑣(𝑟 )𝑖=2𝜋∙𝑟𝑇), para no perder el rango de estabilidad de los
parámetros, de ahí se llega al resto de las soluciones.

> 4.[Resultados:]{.underline}

<table style="width:100%;">
<colgroup>
<col style="width: 7%" />
<col style="width: 7%" />
<col style="width: 7%" />
<col style="width: 7%" />
<col style="width: 7%" />
<col style="width: 7%" />
<col style="width: 7%" />
<col style="width: 7%" />
<col style="width: 7%" />
<col style="width: 7%" />
<col style="width: 7%" />
<col style="width: 7%" />
<col style="width: 7%" />
</colgroup>
<thead>
<tr class="header">
<th colspan="13"><p>El orden de magnitud de las energías y del momento
angular es casi el mismo entre sistema y sistema, siendo el máximo dos
órdenes de diferencia, debido a que los demás parámetros también son de
orden similar.</p>
<p>Se aprecia que la energía potencial es siempre negativa y la cinética
siempre positiva, además</p></th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>de</td>
<td>que</td>
<td>la</td>
<td>energía</td>
<td>total</td>
<td>y</td>
<td>el</td>
<td>momento</td>
<td>angular</td>
<td>total(módulo</td>
<td>y</td>
<td>todas</td>
<td>las</td>
</tr>
<tr class="even">
<td colspan="13">componentes/dirección),calculados con (4) y (5), se
conservan.</td>
</tr>
</tbody>
</table>

1La solución para la que se ha estudiado el error relativo del método
Verlet y Euler es la del sistema Kepler-16.

3

4.1Órbitas coplanarias\
En estos sistemas, solo hay momento angular en la dirección z, por lo
que el momento total en z es el total de todas las direcciones. La razón
de ello es que cuando las órbitas son coplanarias el vector posición y
el de velocidad no tienen en ningún momento componente z (según los ejes
escogidos).

4.1.1 Órbitas externas\
-[Órbitas elípticas, sistema Kepler-16 real]{.underline}\[3\]\[8\]\
La aproximación para la velocidad con la fórmula para órbitas circulares
da buenos resultados, elípticos de todas formas, por lo que parece ser
que un pequeño cambio no afecta al resultado final del equilibrio.

Al ser órbitas periódicas, las energías cinética y potencial en función
del tiempo también lo son. El período de las estrellas es más corto que
el del planeta y cada período de las curvas de energías corresponde a
una órbita entera de las estrellas.

El momento angular total resulta ser positivo, porque los tres cuerpos
giran de forma que el camino más corto desde el vector posición al
vector momento lineal va en dirección contrario a las agujas del reloj.
Es decir, se trata de órbitas 'progrado'\[2\], porque los tres cuerpos
giran en el mismo sentido. Los momentos angulares de cada cuerpo no se
conservan a diferencia del total.

> ![](vertopal_25d73e8d28db483db7232314818a56a9/media/image1.png){width="4.654166666666667in"
> height="3.395832239720035in"}
>
> *Figura 1. Trayectoria, energíasmecánica, cinética y potencialen
> función del tiempo del sistema Kepler-16.*
>
> \-[Órbitas inestables]{.underline}

Al adaptar las condiciones del sistema Kepler para obtener otras órbitas
se dan dificultades para obtener órbitas estables. El planeta debe
estar, para las órbitas externas, a una distancia amplia de las
estrellas, para no adquirir demasiada aceleración. Para las órbitas
internas, sin embargo, la masa debe ser menor y debe estar más cerca de
una de las estrellas, para evitar que la fuerza ejercida por la otra
estrella lo saque de órbita.

4

De no cumplir estas condiciones, el planeta choca con una de las
estrellas o sale despedido debido a una fuerte aceleración, como se
puede observar en la figura 2.

![](vertopal_25d73e8d28db483db7232314818a56a9/media/image2.png){width="6.268054461942257in"
height="3.151388888888889in"}

*Figura 2. Casos no estables. A la izquierda el planeta colapsa contra
una estrella. A la derecha se acerca demasiado a una estrella y sale
despedido.*

> \-[Órbita circular alrededor de estrellas con giros
> internos]{.underline}\[1\]

Modificando las masas, esta vez los tres cuerpos son estrellas2, y,
multiplicando las posiciones y velocidades iniciales por las
proporciones adecuadas\[1\], se llega a estas órbitas.

Las energías cinética y potencial en función del tiempo son
aproximadamente periódicas y las órbitas también. Las energías fluctúan
más que en el caso anterior, debido a que la interacción entre las
estrellas centrales es más fuerte y tardan menos en completar un ciclo
pequeño del aro que las anteriores en completar su órbita.

Como las velocidades de los tres cuerpos son de la misma escala y la
estrella exterior está a gran distancia respecto a las otras, teniendo
en cuenta la ecuación del momento angular es razonable pensar que esta
sea quien más aportar al momento angular total y, efectivamente, el
total es negativo como el de la estrella exterior (a diferencia del de
las otras dos). La de la estrella exterior se trata de una órbita
retrógrada\[2\], que gira en el sentido contrario a la de las estrellas
centrales.

2Nos referimos a estrellas a los cuerpos de masa del orden de la masa
del Sol, y a planetas a cuerpos de masa del orden de la masa de Júpiter.

5

![](vertopal_25d73e8d28db483db7232314818a56a9/media/image3.png){width="7.858332239720035in"
height="3.6236100174978128in"}

> *Figura 3. Trayectoria, energías cinética y potencial (valor absoluto)
> en función del tiempo de la órbita circular alrededor de estrellas con
> giros internos.*
>
> \-[Órbita circular alrededor de dos estrellas con giros internos
> pequeños]{.underline}
>
> Esta figura se consigue igualando la masa de la que era antes más
> pequeña a la de las otras dos. Además, la distancia entre las otras
> dos masas se acorta. Todo ello lleva a dos órbitas de estrellas más
> finas: los dos cuerpos centrales están más juntos, se atraen más, y la
> órbita circular planetaria tiene un radio menor. Las energías cinética
> y potencial en función del tiempo fluctúan aún más en este caso que en
> el anterior, porque las estrellas centrales forman giros internos a su
> órbita principal más pequeños y, por tanto, completan una órbita
> interior en menos tiempo.

![](vertopal_25d73e8d28db483db7232314818a56a9/media/image4.png){width="8.04861111111111in"
height="3.6388877952755907in"}

> *Figura 4. Trayectoria (en 250 días) y energías cinética y potencialen
> función del tiempo de la órbita circular y las dos órbitas estelares
> con giros internos pequeños.*

6

> \-[Órbita planetaria externa]{.underline}\[2\]
>
> Para esta figura se ha minimizado la distancia del planeta a las dos
> estrellas y aumentado la distancia inicial entre las órbitas
> estelares. También, se han utilizado otras proporciones para las
> posiciones y velocidades iniciales\[2\]. De esta forma, las órbitas
> estelares sufren el efecto contrario y casi se solapan. La órbita
> planetaria es más compleja, debido a la cercanía respecto a las
> estrellas.
>
> En este caso, las energías cinética y potencial en función del tiempo
> fluctúan casi con la misma frecuencia de la órbita circular y los dos
> aros con pétalos. Esto es curioso porque, en este caso, la órbita más
> compleja es la del planeta y no la de las estrellas interiores, y a
> pesar de ello, las energías varían de forma parecida, debido a la
> periodicidad y complejidad similar de las órbitas.

![](vertopal_25d73e8d28db483db7232314818a56a9/media/image5.png){width="8.020833333333334in"
height="3.6527777777777777in"}

> *Figura 5. Trayectoria (en 1500 días) y energías cinética y potencial
> (en valor absoluto) en función del tiempo de la órbita planetaria
> externa.*
>
> Como curiosidad, se añadió un segundo centro de masas (de dos cuerpos)
> en varios sistemas para ver qué pasaba y, como era de esperar, este
> orbitaba alrededor del centro de masas común de los tres cuerpos.
>
> 4.1.2Órbitas internas\
> -[Órbita interna alrededor de estrellas con órbita
> circular]{.underline}
>
> Para conseguireste sistema, se parte del código del anterior. Los
> cambios efectuados más notables son las posiciones y velocidades
> iniciales y la masa del planeta. La masa del planeta es 100 veces
> inferior a la de Júpiter, inferior al resto de planetas del trabajo,
> por lo que no afecta demasiado a las órbitas de las estrellas al
> cruzarse con ellas.

7

Las energías cinética y potencial difieren en ∼12.5·1039 J, lo mismo a
efectos prácticos que en los casos anteriores, pero varían mucho menos
en el tiempo, tan poco que parece que las dos se mantienen constantes.
Esto ocurre porque la masa del planeta es muy pequeña, la quinta parte
de la de Júpiter, en comparación con la de las estrellas. Por lo tanto,
las estrellas son las que más energía aportan al sistema y, al ser sus
órbitas muy cercanas a circulares (las estrellas están una en oposición
a la otra, a una distancia aproximadamente constante en el tiempo), sus
energías potenciales y cinética son aproximadamente constantes.

> ![](vertopal_25d73e8d28db483db7232314818a56a9/media/image6.png){width="3.270832239720035in"
> height="2.8402777777777777in"}
>
> *Figura 6. Trayectoria de la órbita planetaria interna alrededor de
> estrellas en órbita circular.*
>
> \-[Órbita interna, planeta orbitando una de las dos
> estrellas]{.underline}\[2\]

Para conseguir este sistema, lo determinante son las posiciones
iniciales, después se ajusta todo para que pueda salir una figura a
partir de esas posiciones. La estrategia consiste en dar un valor
pequeño a la masa, aunque mayor que la real, y ponerla cerca de una de
las estrellas. De esa forma, se asegura que quede confinada al campo
gravitatorio de esa estrella.

Cabe mencionar que, para conseguir el sistema anterior, se empieza
aproximándose por una parecida a esta, ya que las órbitas deseadas
también son internas, y eso otorga los resultados esperados.

En este caso, la variación de las energías cinética y potencial en
función del tiempo se debe a dos razones principales, la órbita del
planeta alrededor de una de las estrellas, y el movimiento relativo de
las estrellas. Como en el caso anterior, las energías las determinan
principalmente las estrellas. En cambio, el planeta está tan cerca de
una de las estrellas que afecta a su energía y crea oscilaciones
pequeñas en las energías totales.

8

![](vertopal_25d73e8d28db483db7232314818a56a9/media/image8.png){width="7.586111111111111in"
height="3.012234251968504in"}

> ![](vertopal_25d73e8d28db483db7232314818a56a9/media/image7.png){width="4.929166666666666in"
> height="3.6486111111111112in"}
>
> *Figura 7. Trayectoria y energías cinética y potencial en valor
> absoluto en función del tiempo del a órbita interna del planeta
> alrededor de una de las dos estrellas.*

+-----------------------------------+-----------------------------------+
| 4.2.                              | > Órbitas no coplanarias\[4\]     |
+===================================+===================================+
+-----------------------------------+-----------------------------------+

> Para obtener un sistema con órbitas no coplanarias a partir del
> sistema Kepler-16, se intercambia la componente y de la estrella
> Kepler-16B con la componente z y, después de ejecutar el código, se
> copian las condiciones del último paso, teniendo en cuenta que los
> sistemas tienden a estabilizarse y estas condiciones son más cercanas
> a las del sistema estable. El momento angular total se distribuye en
> todas sus componentes, pero siendo constante en todas ellas, por lo
> que se conserva en módulo y dirección. La energía mecánica se conserva
> y la cinética y potencial son periódicas en función del tiempo. Aun
> reduciendo la masa del planeta drásticamente, las energías fluctúan
> igual, por lo que la variación de la energía se debe exclusivamente a
> las órbitas elípticas de las estrellas. Cuando las estrellas están más
> cerca entre sí, sufren una aceleración mayor y las energías cinética y
> potencial cambian rápidamente por un tiempo (véanse los mínimos de la
> energía potencial y máximos de la cinética en la figura 8).

9

*Figura 8. Trayectoria y energías mecánica, cinética y potencialen
función del tiempo de lasórbitas no coplanarias.*

> 5.[Conclusiones:]{.underline}

El objetivo de este trabajo era partir del sistema Kepler-16 y llegar a
soluciones estables del problema de los tres cuerpos. Se han logrado
sistemas variados, desde órbitas coplanarias a no coplanarias, tanto
externas como internas y de formas diferentes.

Se ha observado que la periodicidad de las energías cinética y potencial
en función del tiempo está directamente correlacionada con la
periodicidad de las órbitas. Se ha encontrado una relación entre la
forma de las órbitas y las variaciones en las energías estudiadas. Por
otra parte, se han comprobado de forma numérica las leyes de la
conservación tanto de la energía total como la del momento angular total
(módulo y dirección).

Lo normal es que las órbitas de tres cuerpos no sean estables: que no se
atraigan lo suficiente, que un cuerpo sea expulsado y que los otros dos
formen un sistema de dos o que los tres objetos choquen en algún momento
y se alejen infinitamente.

Finalmente, concluir que el problema de los tres cuerpos alberga una
riqueza adicional respecto al de dos, ya que en el de dos las únicas
soluciones estables posibles son las secciones cónicas (órbitas
coplanarias)\[7\]. Se ha estudiado un número limitado de sistemas con la
restricción de partir del sistema Kepler-16 y usando la constante de
gravitación universal. En cambio, la diversidad real del problema es
mucho mayor y todavía es objeto de estudio.

> 6.[Referencias:]{.underline}\
> 1.Sheen, M. (2016, septiembre). *3* *Body* *Solutions*. Github.
>
> 2..
>
> S. Edgeworth.\
> 3.Colaboradores*-16*. Wikipedia, the free Encyclopedia\
> 4.Henríquez, R. A. ital de planetas no coplanares en un sistema
> estelar binario. *Revista de La Escuela de Física*, *5*(1), 1--5.
>
> 5.*od for Integration of Ordinary Differential*
>
> *Equations* *for* *Initial* *Value* *Problems*. Portland State
> University.
>
> 6.niversity of Delaware.
>
> 7.. University of California, Santa Barbara\
> 8.Exoplanam. Kepler-16b. NASA

10

> 7.[Agradecimientos]{.underline}:\
> Querríamos agradecer a la UAM por la licencia de Matlab y a Gabino
> Rubio, nuestro profesor de Computación, por su ayuda en el proyecto,
> tanto a la hora de plantear el código como al buscar información.

11
