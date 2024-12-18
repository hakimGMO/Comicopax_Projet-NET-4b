# Comicopax_Projet-NET-4b
## Description
Projet de couplage de données transcriptomiques, protéomiques et métabolomiques par une approche systémique.
### Contexte Scientifique
- Collaboration : Emmanuelle Becker et Olivier Dameron (Laboratoire INRIA-IRISA) et Florence Gondret (Laboratoire IINRAe-PEGASE)
- Contexte biologique : Intégration de données -omiques de nature différentes
- Méthodes : Parcours de graphes, marches aléatoires
- Format de données : BioPAX

## Fonctionnalités
### Analyse de réseaux protéiques :
- Identification et analyse des chemins les plus courts entre les protéines
- Calcul de centralité de betweenness pour les nœuds du réseau
- Analyse des fréquences d'apparition des nœuds dans les chemins
- Filtrage intelligent des métabolites via une liste noire configurable
- Support des analyses sur graphes dirigés et non-dirigés

### Génération de rapports :
- Rapports HTML interactifs avec tableaux et visualisations
- Rapports texte pour intégration dans des pipelines
- Statistiques détaillées incluant :
  - Vue d'ensemble du réseau
  - Analyse des chemins les plus courts
  - Distribution des fréquences des nœuds
  - Scores de centralité

## Prérequis
- Python 3.7+
- NetworkX
- concurrent.futures (inclus dans Python standard)

## Installation
1. Cloner le dépôt :
```bash
git clone https://github.com/hakimGMO/Comicopax_Projet-NET-4b.git
cd Comicopax_Projet-NET-4b
```

2. Installer les dépendances :
```bash
pip install networkx
```

## Utilisation
Commande de base :
```bash
python protein_network_analyzer.py -f chemin/vers/fichier.graphml
```

Options disponibles :
```
Arguments requis :
  -f, --file           Chemin vers le fichier GraphML

Arguments optionnels :
  -o, --output         Chemin du fichier de sortie (défaut: network_analysis.html)
  --format            Format de sortie ('html' ou 'txt', défaut: html)
  --no-open           Ne pas ouvrir automatiquement le rapport
  -r, --random-nodes  Nombre de nœuds à sélectionner aléatoirement
  -l, --node-list     Liste spécifique de nœuds à analyser
  --disable-blacklist Désactiver la liste noire des métabolites par défaut
```

Exemples d'utilisation :
```bash
# Analyse basique avec rapport HTML
python protein_network_analyzer.py -f network.graphml

# Analyse de nœuds spécifiques avec rapport texte
python protein_network_analyzer.py -f network.graphml -l Protein1 Protein2 --format txt

# Analyse aléatoire avec rapport personnalisé
python protein_network_analyzer.py -f network.graphml -r 10 -o custom_report.html
```

## Structure du projet
```
Comicopax/
├── src/
│   ├── protein_network_analyzer.py
│   ├── utils/
│   └── tests/
├── data/
│   └── example_networks/
├── docs/
├── requirements.txt
└── README.md
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
+ Liens utiles de BIOPAX : https://www.biopax.org
+ Tutoriels sur les bibliothèques de graphes : https://youtube.com/playlist?list=PLGZqdNxqKzfYXTwYAZIlmjnQmrytCSR1J&si=7r3W5AKvKs33zUx0

## Contributeurs
+ Abdelhakim Bouazzaoui
+ Francisco Bretas
+ Corentin Dumortier
+ Pierrick Matelot
+ Raphael Pageau
+ Chenao Wu

## Contact
Pour toute question ou suggestion, veuillez créer une issue sur GitHub ou contacter l'équipe de développement.
