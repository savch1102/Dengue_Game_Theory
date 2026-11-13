# Thresholds and Self-Limiting Cooperation in Awareness-Driven Control of *Aedes*-Borne Disease

A repository for code, scripts and data from the manuscript “Thresholds and self-limiting cooperation in awareness-driven control of *Aedes*-borne disease”, which couples evolutionary game theory with dengue transmission and *Aedes* mosquito dynamics to study how disease awareness and vector awareness shape cooperative breeding-site control.

> Authors: **Santiago A. Villamil-Chacón**, **Jaime Cascante-Vega** & **Mauricio Santos-Vega**.

This work investigates:

+ How awareness of human cases (disease awareness) versus awareness of mosquito abundance (vector awareness) affects the persistence of cooperative vector control.
+ An analytically derived threshold in vector awareness above which cooperation is sustained, and the self-limiting nature of that cooperation (behavior stabilizes well below the level needed to bring $R_0 < 1$).
+ Changes in the dynamical regime (monotonic vs damped-oscillatory return to equilibrium) driven by vector awareness.
+ Interaction of behavioral feedbacks with contrasting seasonal temperature regimes (bimodal Andean vs unimodal Southeast-Asian).

The manuscript contributes to the understanding of when voluntary, awareness-driven prevention can persist between outbreaks, why it tends to plateau short of elimination, and how risk communication focused on vector ecology—and timed to local seasonality—may support community control more effectively than case-based messaging alone. Results are theoretical and should be read as qualitative properties of the coupled system rather than calibrated predictions.

## Repository contents

- Model implementation (ODE system: human SEIR by strategy, mosquito SEI, replicator dynamics).
- Scripts for the (sigma_h, sigma_v) parameter sweeps, eigenvalue analysis, and R0 evaluation.
- Temperature forcing (Anapoima, Colombia; Narathiwat, Thailand), coherence and lag analyses.
- Scripts to reproduce main-text and SI figures.
- Parameter tables and numerical settings (see SI for full derivations).

## Citation

Villamil-Chacón SA, Cascante-Vega J, Santos-Vega M. Thresholds and self-limiting cooperation in awareness-driven control of *Aedes*-borne disease.

## Contact

Corresponding author: Santiago A. Villamil-Chacón — sa.villamil@uniandes.edu.co
