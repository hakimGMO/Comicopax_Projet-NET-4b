# Comicopax_Projet-NET-4b

## Description
Projet de couplage de données transcriptomiques, protéomiques et métabolomiques par une approche systémique.

### Contexte
- Collaboration : Emmanuelle Becker et Olivier Dameron (Laboratoire INRIA-IRISA) et Florence Gondret (Laboratoire IINRAe-PEGASE)
- Contexte biologique : Intégration de données -omiques de nature différentes
- Méthodes : Parcours de graphes, marches aléatoires
- Format : BioPAX

## Fonctionnalités
- Analyse de réseaux protéiques
- Calcul des chemins les plus courts
- Analyse de centralité des nœuds
- Identification des nœuds ubiquitaires
- Génération de rapports HTML et TXT
- Support du traitement parallèle

## Prérequis
- Python 3.x
- NetworkX

## Installation
1. Cloner le dépôt :
```bash
git clone https://github.com/hakimGMO/Comicopax_Projet-NET-4b.git
cd Comicopax_Projet-NET-4b
```

2. installer les dépendances
```bash
pip install networkx
```
## Utilisation
Execution de base :
```bash
python protein_network_analyzer.py -f chemin/vers/fichier.graphml
```
Options possibles :
Le chemin vers le fichier GraphML qui est obligatoire.
```bash
-f ou --file
```
Le nombre d'itérations pour l'analyse des chemins par défaut mis à 100 itérations.
```bash
-i ou --iterations 
```
Le nombre d'échantillons pour le calcul de centralité par défaut mis à 100 échantillons.
```bash
-k ou --centrality-samples
```
Le chemin du fichier de sortie qui sera par défaut "./network_analysis.html".
```bash
-o ou --output
```
Le format de sortie du rapport HTML ou TXT, qui est par défaut un rapport HTML.
```bash
--format
```
Permet de ne pas ouvrir automatiquement le rapport à la fin de l'éxecution.
```bash
--no-open
```
le seuil pour l'identification des nœuds ubiquitaires par défaut mis à 0.9.
```bash
-t ou --threshold
```
## Workflow git
Mettre a jour son dépôt local :
```bash
git pull
```
Soumettre des modifications :
```bash
git commit -a -m "Description des modifications"
git push
```
## Documentation
+ Documentation NetworkX : https://networkx.org/documentation/stable/reference/readwrite/generated/networkx.readwrite.graphml.read_graphml.html
+ Tutoriels sur les bibliothèques de graphes : https://youtube.com/playlist?list=PLGZqdNxqKzfYXTwYAZIlmjnQmrytCSR1J&si=7r3W5AKvKs33zUx0

## Contributeurs

+ Abdelhakim Bouazzaoui
+ Francisco Bretas
+ Corentin Dumortier
+ Pierrick Matelot
+ Raphael Pageau
+ Chenao Wu
