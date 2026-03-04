# MPNA Project – TeaLeaf Benchmark Analysis

## Description

Ce dépôt contient les résultats du projet **MPNA – Analyse du benchmark TeaLeaf**.

L'objectif du projet est :

- étudier le fonctionnement du mini-app **TeaLeaf**
- analyser ses performances CPU / GPU
- implémenter un **mini solveur inspiré de TeaLeaf**
- réaliser une analyse de performance avec **MAQAO**

Le code fourni (`miniTealeaf2D`) est une implémentation simplifiée en **Fortran** d’un solveur de diffusion sur **grille 2D** utilisant la méthode **Gradient Conjugate (CG)**.

---

# Structure du projet

```
reportMPNA-Tealeaf/
│
├── RapportMPNA_Kim_Selyan_Franc.pdf
│   Rapport complet du projet
│
├── PresentationMPNA_Kim_Selyan_Franc.pdf
│   Support de présentation
│
├── miniTealeaf2D/
│   Code source Fortran du mini solveur
│
│   ├── main.f90
│   ├── config.f90
│   ├── hydro.f90
│   ├── radiation.f90
│   ├── io.f90
│   ├── Makefile
│   ├── miniTealeaf2D
│   └── plot
│
└── analyseBenchMakao/
    Analyse de performance avec MAQAO
    │
    ├── analyseMaqao.ipynb
    ├── analyseMaqao.pdf
    ├── reportsMaqao.csv
    └── figures/
```

---

# Compilation

Se placer dans le dossier du code :

```bash
cd miniTealeaf2D
```

Compiler :

```bash
make
```

Cela génère l’exécutable :

```
miniTealeaf2D
```

---

# Exécution

Lancer la simulation :

```bash
./miniTealeaf2D
```

Le programme :

1. initialise une grille 2D  
2. initialise un champ d’énergie  
3. résout l’équation de diffusion avec un solveur **Gradient Conjugate**  
4. écrit les résultats dans les fichiers de sortie  

---

# Visualisation (optionnel)

Un script de visualisation est fourni.

Depuis le dossier `miniTealeaf2D` :

```bash
./plot
```

Le script utilise **gnuplot** pour afficher l’évolution des données produites par la simulation.

---

# Analyse de performance

Le dossier :

```
analyseBenchMakao/
```

contient l’analyse réalisée avec **MAQAO** :

- notebook d’analyse
- rapport PDF
- données CSV
- figures générées

---

# Environnement de test

Le code a été testé avec :

```
Linux
gfortran
gnuplot
```

---

# Auteurs

Projet réalisé par :

Kim Davy NDAYONGEJE  
Selyan KABLIA  
Franc Chanel TCHOLIEU NSAMOU  

Master CHPS : UVSQ/ENS – Université Paris-Saclay