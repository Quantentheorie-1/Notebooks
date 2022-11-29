# Quantentheorie I

Dieses Repository enthält Jupyter Notebooks, die zur Begleitung der
**Quantentheorie I Vorlesung mit Übung** an der Technischen Universität Wien zur Verfügung gestellt
werden. Die frei verfügbaren Jupyter Notebooks dienen der Visualisierung und Durchführung
von "numerischen Experimenten", wodurch ein besseres Verständnis
der Quantentheorie ermöglicht werden soll. Es wird vorausgesetzt das Studierende Parallel dazu,
die Vorlesung und Rechenübungen besuchen.

Die verwendeten Algorithmen und numerischen Verfahren erlauben es auch analytisch nicht exakt lösbare
Fragestellungen zu komplexen quantenmechanischen Systemen zu untersuchen.
Studierende werden ausdrücklich dazu ermutigt die verfügbaren Notebooks und gewählten Parameter
zu modifizieren um diverse Fragestellungen mittels "numerischer Experimente" zu untersuchen. 

## Dynamik eindimensionaler Wellenpakete

In diesem Kapitel betrachten wir die zeitliche Änderung eindimensionaler Wellenpakete,
welche sich aus der Lösung der zeitabhängigen Schrödingergleichung ergibt.
In den folgenden Abschnitten wird die zeitabhängige Schrödingergleichung,

<img src="https://render.githubusercontent.com/render/math?math=i\hbar \frac{\partial}{\partial t} \Psi(x,t)=\hat{H} \Psi(x,t)">,

für unterschiedliche Fälle numerisch gelöst. Eine eindeutige Lösung der obigen Differentialgleichung für
<img src="https://render.githubusercontent.com/render/math?math=t \ge t_0">
setzt voraus, dass
<img src="https://render.githubusercontent.com/render/math?math=\hat{H}"> und
<img src="https://render.githubusercontent.com/render/math?math=\Psi(x,t_0)"> gegeben ist.

Die zeitabhängige Wellenfunktion für *zeitunabhängige* Hamiltonoperatoren folgt aus
folgender unitärer Transformation

<img src="https://render.githubusercontent.com/render/math?math=\Psi(x,t)= e^{-\frac{i}{\hbar}\hat{H}(t-t_0)} \Psi(x,t_0)">.

Diese Transformation wird auch Zeitentwicklung oder Zeitpropagation genannt.
Es existieren viele unterschiedliche numerische Verfahren zur Berechnung
der zeitabhängigen Wellenfunktion.

### Split-Operator-Methode

In den folgenden Abschnitten verwenden wir die
[Split-Operator-Methode](https://de.wikipedia.org/wiki/Split-Operator-Methode) um die Wirkung des
Zeitentwicklungsoperators <img src="https://render.githubusercontent.com/render/math?math=e^{-\frac{i}{\hbar} \hat{H}(t-t_0)}">
auf die Wellenfunktion zu vereinfachen und näherungsweise zu berechnen.

Eine wichtige Voraussetzung der Split-Operator-Methode ist, dass der Hamiltonoperator des
betrachteten Systems durch die Summe eines kinetischen Energieoperators
<img src="https://render.githubusercontent.com/render/math?math=\hat{T}=-\frac{\hbar^2}{2 m}\Delta">
und eines lokalen Potentialoperators
<img src="https://render.githubusercontent.com/render/math?math=\hat{V}=V(x)">
gegeben ist.
Damit ergibt sich folgende zentrale Näherung der Split-Operator-Methode für den Zeitentwicklungsoperator

<img src="https://render.githubusercontent.com/render/math?math=e^{-\frac{i}{\hbar} \hat{H}\Delta t} \approx e^{-\frac{i}{\hbar} \frac{\hat{T}}{2}\Delta t}  e^{-\frac{i}{\hbar} \hat{V} \Delta t} e^{-\frac{i}{\hbar} \frac{\hat{T}}{2}\Delta t} ">,

welche einen numerischen Fehler mit
<img src="https://render.githubusercontent.com/render/math?math=\mathcal{O}(\Delta t)^3">
verursacht. Der numerische Fehler folgt aus dem Vergleich der Reihenentwicklung nach
<img src="https://render.githubusercontent.com/render/math?math=\Delta t">
von
<img src="https://render.githubusercontent.com/render/math?math=e^{-\frac{i}{\hbar} \hat{H}\Delta t}">
und 
<img src="https://render.githubusercontent.com/render/math?math=e^{-\frac{i}{\hbar} \frac{\hat{T}}{2}\Delta t}  e^{-\frac{i}{\hbar} \hat{V} \Delta t} e^{-\frac{i}{\hbar} \frac{\hat{T}}{2}\Delta t} ">.
Hinreichend kleine Zeitschritte führen zu vernachlässigbar kleinen numerischen
Fehlern in  <img src="https://render.githubusercontent.com/render/math?math=\Psi(x,t)">.

In der Split-Operator-Methode wird die Wirkung der beiden unterschiedliechen Terme des
genäherten Zeitentwicklungsoperators getrennt und in unterschiedlichen Darstellungen
berechnet.
Die Wirkung von 
<img src="https://render.githubusercontent.com/render/math?math=e^{-\frac{i}{\hbar}\frac{\hat{T}}{2}\Delta t}">
wird im Impulsraum berechnet, während die Wirkung von
<img src="https://render.githubusercontent.com/render/math?math=e^{-\frac{i}{\hbar} \hat{V} \Delta t}">
im Ortsraum berechnet wird.

### Freies Teilchen (vgl. 3.2.1 im Skriptum)

Wir betrachten zuerst die Zeitentwicklung eines Gaußschen Wellenpaketes für den Fall
eines freien Teilchens (*V=0*).
In diesem Fall ist der Zeitentwicklungsoperator gegeben durch:
<img src="https://render.githubusercontent.com/render/math?math=e^{-\frac{i}{\hbar}\hat{T}(t-t_0)}">,
welcher sich im Impulsraum zu
<img src="https://render.githubusercontent.com/render/math?math=e^{-\frac{i}{\hbar}\frac{\hbar^2 k^2}{2 m}\Delta t}">
vereinfacht.

Wir definieren das folgende Gaußsche Wellenpaket im Impulsraum zum Zeitpunlt *t=0*:
<img src="https://render.githubusercontent.com/render/math?math=\tilde{\Psi}(k)= e^{-(k-k_0)^2/d^2}">.
Der mittlere Impuls dieses Wellenpaketes ist gegeben durch
<img src="https://render.githubusercontent.com/render/math?math=\hbar k_0">.

Das Jupyter Notebook im untenstehenden Link löst die zeitabhänige Schrödingergleichung für das freie
Teilchen und visualisiert das zeitabhängige Verhalten
- der Wellenfunktion im Impulsraum <img src="https://render.githubusercontent.com/render/math?math=\tilde{\Psi}(k,t)">.
- der Wellenfunktion im Ortsraum <img src="https://render.githubusercontent.com/render/math?math={\Psi}(r,t)">.
- der Aufenthaltswahrscheinlichkeitsdichte im Ortsraum <img src="https://render.githubusercontent.com/render/math?math=|{\Psi}(r,t)|^2">.

Es ist zu beobachten, dass
- sich der Schwerpunkt von <img src="https://render.githubusercontent.com/render/math?math=|{\Psi}(r,t)|^2"> nach rechts bewegt. 
- <img src="https://render.githubusercontent.com/render/math?math=|{\Psi}(r,t)|^2"> mit der Zeit "zerfließt".
- <img src="https://render.githubusercontent.com/render/math?math={\tilde{\Psi}}(k,t)"> bei größeren Werten von *k* schneller "oszilliert".

Versuchen sie nun folgende Fragen mit analytischen Überlegungen und numerischen Experimenten zu beantworten:
- Was beobachten sie, wenn sie den Parameter *d* vergrößern/verkleinern?
- Welche Aussagen können Sie über mögliche Messwerte des Aufenthaltsorts des Teilchens als Funktion der Zeit treffen?
- Welche Aussagen können Sie über mögliche Messwerte des Impulses des Teilchenis als Funktion der Zeit treffen?
- Was passiert im Grenzfall <img src="https://render.githubusercontent.com/render/math?math=d \rightarrow 0"> ?
- Was passiert im Grenzfall <img src="https://render.githubusercontent.com/render/math?math=d \rightarrow \infty"> ?

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Quantentheorie-1/Notebooks/blob/main/notebooks/TD-Free-Schroedinger.ipynb)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Quantentheorie-1/Notebooks/HEAD?labpath=notebooks%2FTD-Free-Schroedinger.ipynb)



### Teilchen an der Potentialbarriere (vgl. 3.4 im Skriptum)


[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Quantentheorie-1/Notebooks/blob/main/notebooks/TD-Schroedinger-Barriere.ipynb)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Quantentheorie-1/Notebooks/HEAD?labpath=notebooks%2FTD-Schroedinger-Barriere.ipynb)


### Leiteroperatoren zur algebraischen Lösung des harmonischen Oszillators (vgl. 5.2.1 im Skriptum)

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Quantentheorie-1/Notebooks/blob/main/notebooks/Ladder-Operators.ipynb)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Quantentheorie-1/Notebooks/HEAD?labpath=notebooks%2FLadder-Operators.ipynb)


### Kohärente Glauber-Zustände (vgl. 5.3.2 im Skriptum)

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Quantentheorie-1/Notebooks/blob/main/notebooks/TD-Schroedinger-Glauber.ipynb)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Quantentheorie-1/Notebooks/HEAD?labpath=notebooks%2FTD-Schroedinger-Glauber.ipynb)



