# Quantentheorie I

Dieses Repository enthält Jupyter Notebooks, die zur Begleitung der
**Quantentheorie I Vorlesung mit Übung** an der Technischen Universität Wien zur Verfügung gestellt
werden. Die frei verfügbaren Jupyter Notebooks dienen der Visualisierung und Durchführung
von "numerischen Experimenten", wodurch ein besseres Verständnis
der Quantentheorie ermöglicht werden soll. Es wird vorausgesetzt das Studierende Parallel dazu
die Vorlesung und Rechenübungen besuchen.

Die verwendeten Algorithmen und numerischen Verfahren erlauben es auch analytisch nicht exakt lösbare
Fragestellungen zu komplexen quantenmechanischen Systemen zu untersuchen.
Studierende werden ausdrücklich dazu ermutigt die verfügbaren Notebooks und gewählten Parameter
zu modifizieren um diverse Fragestellungen mittels "numerischer Experimente" zu untersuchen. 

## Dynamik eindimensionaler Wellenpakete

In diesem Kapitel betrachten wir die zeitliche Änderung eindimensionaler Wellenpakete,
welche sich aus der Lösung der zeitabhängigen Schrödingergleichung ergibt.
In den folgenden Abschnitten wird die zeitabhängige Schrödingergleichung,

$$ i\hbar \frac{\partial}{\partial t} \Psi(x,t)=\hat{H} \Psi(x,t), $$

für unterschiedliche Fälle numerisch gelöst. Eine eindeutige Lösung der obigen Differentialgleichung für 
$t \ge t_0$ 
setzt voraus, dass
$\hat{H}$ und
$\Psi(x,t_0)$ gegeben ist.

Die zeitabhängige Wellenfunktion für *zeitunabhängige* Hamiltonoperatoren folgt aus
folgender unitärer Transformation

$$\Psi(x,t)= e^{-\frac{i}{\hbar}\hat{H}(t-t_0)} \Psi(x,t_0).$$

Diese Transformation wird auch Zeitentwicklung oder Zeitpropagation genannt.
Es existieren viele unterschiedliche numerische Verfahren zur Berechnung
der zeitabhängigen Wellenfunktion.

### Split-Operator-Methode

In den folgenden Abschnitten verwenden wir die
[Split-Operator-Methode](https://de.wikipedia.org/wiki/Split-Operator-Methode) um die Wirkung des
Zeitentwicklungsoperators $e^{-\frac{i}{\hbar} \hat{H}(t-t_0)}$
auf die Wellenfunktion zu vereinfachen und näherungsweise zu berechnen.

Eine wichtige Voraussetzung der Split-Operator-Methode ist, dass der Hamiltonoperator des
betrachteten Systems durch die Summe eines kinetischen Energieoperators
$\hat{T}=-\frac{\hbar^2}{2 m}\Delta$
und eines lokalen Potentialoperators
$\hat{V}=V(x)$
gegeben ist.
Damit ergibt sich folgende zentrale Näherung der Split-Operator-Methode für den Zeitentwicklungsoperator

$$e^{-\frac{i}{\hbar} \hat{H}\Delta t} \approx e^{-\frac{i}{\hbar} \frac{\hat{T}}{2}\Delta t}  e^{-\frac{i}{\hbar} \hat{V} \Delta t} e^{-\frac{i}{\hbar} \frac{\hat{T}}{2}\Delta t} ,$$

welche einen numerischen Fehler mit
$\mathcal{O}(\Delta t)^3$
verursacht. Der numerische Fehler folgt aus dem Vergleich der Reihenentwicklung nach
$\Delta t$
von
$e^{-\frac{i}{\hbar} \hat{H}\Delta t}$
und 
$e^{-\frac{i}{\hbar} \frac{\hat{T}}{2}\Delta t}  e^{-\frac{i}{\hbar} \hat{V} \Delta t} e^{-\frac{i}{\hbar} \frac{\hat{T}}{2}\Delta t} .$
Hinreichend kleine Zeitschritte führen zu vernachlässigbar kleinen numerischen
Fehlern in  $\Psi(x,t)$ .

In der Split-Operator-Methode wird die Wirkung der beiden unterschiedliechen Terme des
genäherten Zeitentwicklungsoperators getrennt und in unterschiedlichen Darstellungen
berechnet.
Die Wirkung von 
$e^{-\frac{i}{\hbar}\frac{\hat{T}}{2}\Delta t}$
wird im Impulsraum berechnet, während die Wirkung von
$e^{-\frac{i}{\hbar} \hat{V} \Delta t}$
im Ortsraum berechnet wird.

### Freies Teilchen (vgl. 3.2.1 im Skriptum)

Wir betrachten zuerst die Zeitentwicklung eines Gaußschen Wellenpaketes für den Fall
eines freien Teilchens ($V=0$).
In diesem Fall ist der Zeitentwicklungsoperator gegeben durch:
$e^{-\frac{i}{\hbar}\hat{T}(t-t_0)}$,
welcher sich im Impulsraum zu
$e^{-\frac{i}{\hbar}\frac{\hbar^2 k^2}{2 m}\Delta t}$
vereinfacht.

Wir definieren das folgende Gaußsche Wellenpaket im Impulsraum zum Zeitpunlt $t=0$:
$\tilde{\Psi}(k)= e^{-(k-k_0)^2/d^2}$.
Der mittlere Impuls dieses Wellenpaketes ist gegeben durch
$\hbar k_0$.

Das Jupyter Notebook im untenstehenden Link löst die zeitabhänige Schrödingergleichung für das freie
Teilchen und visualisiert das zeitabhängige Verhalten
- der Wellenfunktion im Impulsraum $\tilde{\Psi}(k,t)$.
- der Wellenfunktion im Ortsraum ${\Psi}(r,t)$.
- der Aufenthaltswahrscheinlichkeitsdichte im Ortsraum $|{\Psi}(r,t)|^2$.

Es ist zu beobachten, dass
- sich der Schwerpunkt von $|{\Psi}(r,t)|^2$ nach rechts bewegt. 
- $|{\Psi}(r,t)|^2$ mit der Zeit "zerfließt".
- ${\tilde{\Psi}}(k,t)$ bei größeren Werten von $k$ schneller "oszilliert".

Versuchen sie nun folgende Fragen mit analytischen Überlegungen und numerischen Experimenten zu beantworten:
- Was beobachten sie, wenn sie den Parameter $d$ vergrößern/verkleinern?
- Welche Aussagen können Sie über mögliche Messwerte des Aufenthaltsorts des Teilchens als Funktion der Zeit treffen?
- Welche Aussagen können Sie über mögliche Messwerte des Impulses des Teilchenis als Funktion der Zeit treffen?
- Was passiert im Grenzfall $d \rightarrow 0$ ?
- Was passiert im Grenzfall $d \rightarrow \infty$ ?

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Quantentheorie-1/Notebooks/blob/Projektarbeit_Ranner/notebooks/TD-Free-Schroedinger.ipynb)


### Teilchen an der Potentialbarriere (vgl. 3.4 im Skriptum)

Nun betrachten wir die Zeitentwicklung eines Gaußschen Wellenpaketes, welches auf eine Potentialbarriere trifft. 
Die Potentialbarriere definieren wir als Gaußsche Glockenkurve:

$$V(x)= v_0  e^{-(x-x_0)^2/b^2}$$

mit $v_0$ als Höhe der Potentialbarriere, $x_0$ als die Position des Maximums des Potentials und $b$ als die Breite der Potentialbarriere. Damit ist auch der lokale Potentialoperator gegeben als

$$\hat{V} = v_0  e^{-(x-x_0)^2/b^2}.$$

Die Wellenfunktion wird analog zum Fall des freien Teilchens initialisiert. 

Der Zeitentwicklungsoperator $e^{-\frac{i}{\hbar}\hat{T}(t-t_0)}$ wird nun wie oben beschrieben aufgeteilt. 

$$e^{-\frac{i}{\hbar} \hat{H}\Delta t} \approx e^{-\frac{i}{\hbar} \frac{\hat{T}}{2}\Delta t}  e^{-\frac{i}{\hbar} \hat{V} \Delta t} e^{-\frac{i}{\hbar} \frac{\hat{T}}{2}\Delta t}$$

Der Potentialterm des Zeitentwicklungsoperators nimmt nun im Ortsraum die Form 
$e^{-\frac{i}{\hbar} \hat{V}} = e^{-\frac{i}{\hbar} v_0 e^{-(x-x_0)^2/b^2}}$ an. 
Die zwei kinetischen Terme sind im Impulsraum gegeben als 
$e^{-\frac{i}{\hbar}\frac{\hat{T}}{2}\Delta t} = e^{-\frac{i}{\hbar}\frac{\hbar^2 k^2}{2 m}\frac{\Delta t}{2}}$.

Damit diese Terme auf die Wellenfunktion in der jeweils passenden Darstellung wirken können, muss auch die Wellenfunktion im entsprechenden Raum gegeben sein. Daher muss die Wellenfunktion mithilfe von Fourier-Transformationen in die jeweilige Darstellung gebracht werden. Insgesamt stellt sich die Anwednung das Zeitentwicklungsoperators auf $\tilde{\Psi}(k,t)$ also folgendermaßen dar:
- Anwendung von $e^{-\frac{i}{\hbar}\frac{\hat{T}}{2}\Delta t}$
- Inverse Fourier Transformation in den Ortsraum
- Anwendung von $e^{-\frac{i}{\hbar} \hat{V}}$
- Fourier Transformation in den Impulsraum
- Anwendung von $e^{-\frac{i}{\hbar}\frac{\hat{T}}{2}\Delta t}$

Es ist zu beobachten, dass
- ein Teil der Wellenfunktion reflekiert und ein Teil transmitiert wird.
- es im Bereich des Potentials zu Interferenzen zwischen den Anteilen mit unterschiedlicher Bewegungsrichtung kommt.

Versuchen sie nun folgende Fragen mit analytischen Überlegungen und numerischen Experimenten zu beantworten:
- Was beobachten sie, wenn sie den Parameter $v_0$ vergrößern/verkleinern?
- Was passiert im Grenzfall $v_0 \rightarrow 0$ ?
- Was passiert im Grenzfall $v_0 \rightarrow \infty$ ?
- Vergleichen sie die Wellenfunktion im Impulsraum vor und nach der Streuung. Was können sie daraus für die Bewegung der Wellenfunktion im Ortsraum schließen?

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Quantentheorie-1/Notebooks/blob/Projektarbeit_Ranner/notebooks/TD-Schroedinger-Barriere.ipynb)

### Kohärente Glauber-Zustände (vgl. 5.3.2 im Skriptum)

Hier betrachten wir die zeitliche Entwicklung von kohärenten Glauber-Zuständen im Potential $V(x)=\frac{x^2}{2}$. Für $t=0$ ist die Wellenfunktion des kohärenten Zustandes im Ortsraum gegeben als:

$$
\Psi(x,0) = \frac{1}{2} e^{ipx/\hbar}e^{-(x)^2/(2\hbar)}
$$

$p$ bezeichnet hier den Erwartungswert des Impulses.

Kohärente Glauber-Zustände sind als Eigenzustände des Absteigeoperators definiert. Der Absteigeoperator im Ortsraum ist gegeben als

$$
\hat{a} = \frac{1}{\sqrt{2}} \left( \sqrt{\frac{m\omega}{\hbar}} x + \sqrt{\frac{\hbar}{m\omega}} \frac{d}{dx}\right)
$$

Lässt man diesen nun auf $\Psi(x,0)$ wirken, erhält man:

$$
\hat{a}\Psi(x,0) = \frac{i p}{\sqrt{2}\hbar} \frac{1}{2} e^{ipx/\hbar}e^{-(x)^2/(2\hbar)} 
= \frac{i p}{\sqrt{2}\hbar} \Psi(x,0)
$$

Somit ist gezeigt, dass $\Psi(x,0)$ tatsächlich ein Glauber-Zustand ist. 

Danach wird analog zum Fall des Teilchens an der Potentialbarriere die Zeitentwicklung des Zustandes im Potential berechnet. 

Es ist zu beobachten, dass
- die Wellenfunktion nicht "zerfließt".
- das Maximum der Aufenthaltswahrscheinlichkeit sich bewegt wie im Fall des klassischen harmonischen Oszillators.

Was beobachten sie, wenn sie den Parameter $p$ vergrößern/verkleinern?
- Was passiert im Grenzfall $p \rightarrow 0$ ?
- Was passiert im Grenzfall $p \rightarrow \infty$ ? Denken sie daran, dass es sich hier um eine numerische Simulation handelt.

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Quantentheorie-1/Notebooks/blob/Projektarbeit_Ranner/notebooks/TD-Schroedinger-Glauber.ipynb)


## Leiteroperatoren zur algebraischen Lösung des harmonischen Oszillators (vgl. 5.2.1 im Skriptum)

Wir visualisieren nun die Wirkung der Leiteroperatoren auf die Wellenfunktionen der Lösungen des harmonischen Oszillators mit Potential $V(x) = \frac{1}{2}m\omega^2 x^2$. Dafür starten wir mit der Wellenfunktion im Grundzustand $\psi_0 (x)$, welche (auf eins normiert) folgendermaßen gegeben ist:

$$
\psi_0 (x) = \frac{m\omega}{\sqrt{2\pi\hbar^2} } e^{-\frac{1}{2} \frac{m\omega}{\hbar} x^2}
$$

Auf diese Wellenfunktion können wir nun die Aufsteige- und Absteigeoperatoren wirken lassen, welche im Ortsraum folgendermaßen definiert sind:

$$
\hat{a}^{\dagger} = \frac{1}{\sqrt{2}} \left( \sqrt{\frac{m\omega}{\hbar}} x - \sqrt{\frac{\hbar}{m\omega}} \frac{d}{dx}\right)
$$

$$
\hat{a} = \frac{1}{\sqrt{2}} \left( \sqrt{\frac{m\omega}{\hbar}} x + \sqrt{\frac{\hbar}{m\omega}} \frac{d}{dx}\right)
$$

Dargestellt wird jeweils die Wellenfunktion, sowie die Aufenthaltswahrscheinlichkeit $\rho (x) = |\psi(x)|^2$.

Versuchen sie nun folgende Fragen mit analytischen Überlegungen und numerischen Experimenten zu beantworten:
- Prüfen sie anhand der Graphen, ob $\hat{a}$ und $\hat{a}^{\dagger}$ kommutieren.
- Was passiert, wenn sie den Absteigeoperator auf den Grundzustand anwenden?
- Ist es möglich, nach Anwendung von $\hat{a}$ auf den Grundzustand, diesen durch Anwendung von $\hat{a}^{\dagger}$ wiederherzustellen?

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Quantentheorie-1/Notebooks/blob/Projektarbeit_Ranner/notebooks/Ladder-Operators.ipynb)


## Kugelflächenfunktionen (vgl 6.3.2 im Skriptum)

Mit diesem Notebook können sie die Kugelflächenfunktionen für $l \in [0,3]$ als interaktives 3D plot darstellen. Die Kugelflächenfunktion sind durch 

$$
Y_l^m(\theta, \phi) = \frac{1}{\sqrt{2\pi}}N_l^m P_l^m\left(cos\theta\right) e^{im\phi}
$$

gegeben, mit

$$
N_l^m = \sqrt{\frac{2l+1}{2}\frac{(l-m)!}{(l+m)!}},
$$

$$
P_l^m(x) = \frac{(-1)^m}{2^l l!} (1-x^2)^{\frac{m}{2}} \frac{d^{l+m}}{dx^{l+m}}(x^2-1)^l.
$$

Gezeigt werden zwei unterschiedliche Darstellungsformen des Realteiles von $Y_l^m$. In der ersten Darstellung gibt die Entfernung eines Flächenpunktes vom Ursprung den Betrag von $Re(Y_l^m)$ in die jeweilige Raumrichtung an. Raumrichtungen, in denen $Re(Y_l^m)$ positiv ist, werden mit rot gekennzeichnet, und jene wo $Re(Y_l^m)$ negativ ist mit blau. 

In der zweiten Darstellung wird der Wert von $Re(Y_l^m)$ mithilfe einer colormap auf der Einheitskugel dargestellt. Nulldurchgängen werden durch eine schwarzen Linie hervorgehoben.

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Quantentheorie-1/Notebooks/blob/Projektarbeit_Ranner/notebooks/Kugelflächenfunktionen.ipynb)

# Quantentheorie II

## Harmonische Störung (vgl. 5.4.4 im Skriptum)

Wir betrachten eine harmonische Störung in zeitabhängiger Störungstheorie. Zum Zeitpunkt $t=0$ wird ein harmonisches Störpotential der Form $\hat{V}(t)=\hat{V} _ {0} sin(\omega t)$ eingeschalten. Für den Übergang vom Zustand $|\phi_i\rangle$ in den Zustand $|\phi_f\rangle$ (mit Energiedifferenz $\epsilon_{if}=\hbar\omega_{if}$) ist die Übergangswahrscheinlichkeit in erster Ordnung, $P^{(1)}_ {if}$, folgendermaßen gegeben:

$$
P^{(1)}_ {if} \propto 
\frac{sin^2(\frac{\omega_{if}+\omega}{2}t)}{\frac{(\omega_{if}+\omega)^2}{4}} + 
\frac{sin^2(\frac{\omega_{if}-\omega}{2}t)}{\frac{(\omega_{if}-\omega)^2}{4}} - 
2 \frac{sin(\frac{\omega_{if}+\omega}{2}t) sin(\frac{\omega_{if}-\omega}{2}t)}
{\frac{(\omega_{if}+\omega)(\omega_{if}-\omega)}{4}} cos(2\omega t)
$$

Die Übergangswahrscheinlichkeit ist von drei Variablen ($t$, $\omega$ und $\omega_{if}$) abhängig. Zur Visualisierung wird die zeitliche Abhängigkeit von $P_{if}^{(1)}(\omega)$ für festgehaltenes $\omega_{if}$ und von $P_{if}^{(1)}(\omega_{if})$ für festgehaltenes $\omega$ animiert.

- Vergleichen sie bei $P_{if}^{(1)}(\omega_{if})$ die Interferenzen zwischen den beiden Peaks für verschiedene Werte von $\omega$.
- Man beachte die lange Periodizität in $t$, $P_{if}^{(1)}(\omega_{if})$ wird erst ungefähr bei $t\approx 1000$ wieder konstant null, bei $P_{if}^{(1)}(\omega)$ passiert dies erst nach doppelt so langer Zeit.

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Quantentheorie-1/Notebooks/blob/Projektarbeit_Ranner/notebooks/Harmonische_Störung.ipynb)
