# Winning project of the Respira Hackathon

# Cross-Cohort Lung Phenotypes after Severe Respiratory Infection

**Digital phenotyping of post-infectious sequelae (Challenge 2).** We looked for groups of patients with lasting sequelae after a severe respiratory infection (mostly COVID-19 ICU and ward patients) that **hold up when you move to a different cohort**, described **how they evolve over time**, and studied **when in follow-up they can be predicted**.

The project works with a **dataset of real-world clinical data from real patients**, provided under privacy safeguards during the hackathon. It merges four Spanish clinical cohorts (CIBERESUCICOVID, POSTCOVID-Lleida, TENACITY and Virgen del Rocío) and covers **9,809 unique patients**: demographics, comorbidities, hospital and ICU stay, and follow-up visits with lung function tests, symptoms, CT imaging and questionnaires.

> **The question was not "are there subgroups?"** Clustering always finds subgroups. The real question is whether they **survive a change of cohort**. The whole project is built around that idea.

---

## Key results

### 1. A phenotype that travels across cohorts
Using only what all cohorts measure the same way (**DLCO, FVC and FEV1** at ~3 months after discharge), two phenotypes emerge:

| Phenotype | Patients | DLCO | FVC | FEV1 |
|---|---|---|---|---|
| **F1 · Preserved function** | 821 | 79 % | 96 % | 99 % |
| **F2 · Functional impairment** | 753 | 65 % | 73 % | 76 % |

- **Stable:** bootstrap Jaccard 0.92.
- **Not a cohort artefact:** Cramér's V between phenotype and cohort is 0.12.
- **Replicates out of sample:** trained on two cohorts and tested on the third, ARI 0.50–0.73.
- **Negative result, reported anyway:** richer symptom-based phenotypes found in a single cohort (CIBERESUCICOVID) **did not replicate** elsewhere.

### 2. The phenotypes evolve differently
Both groups improve, but the gap barely closes. A linear mixed model of DLCO over time gives:

| | 3 months | 6 months | 12 months |
|---|---|---|---|
| F1 · Preserved function | 79.4 % | 82.7 % | 86.1 % |
| F2 · Functional impairment | 65.4 % | 70.0 % | 74.6 % |
| **Difference F2 − F1** | −14.0 | −12.7 | **−11.5** [−13.6, −9.3] |

At one year F2 is still below the 80 % clinical threshold on average. The gap replicates in each cohort and holds after correcting for informative dropout (IPW) and for regression to the mean. In POSTCOVID-Lleida it persists for up to 4 years.

### 3. Prediction works at the 3-month visit, not at discharge
- **At discharge**, every model tops out at **AUC ≈ 0.63** (logistic regression and EBM give the same result). Adding acute-phase ICU data barely helps (≈ 0.65). The limit is the information available at discharge, not the algorithm.
- **At the 3-month visit**, **DLCO alone** predicts who will still be impaired at one year: **AUC 0.79** among impaired patients in unseen cohorts. **No more complex model beats it.**

**Pocket table:** of the patients with DLCO < 80 % at 3 months, how many recover by one year?

| DLCO at 3 months | Recover at 1 year |
|---|---|
| < 60 % | 7 % [4–13] |
| 60–69 % | 24 % [17–33] |
| 70–79 % | 62 % [53–71] |

**Clinical takeaway:** the spirometry + DLCO test at 3 months should guide follow-up. F1 tends to normalise by around 6 months. F2 needs follow-up beyond the first year.

---

## Approach: the phenotype "passport"

Before running any clustering we fixed, in [`config.yaml`](config.yaml), what a phenotype must pass to count as real. A phenotype that fails is still reported, as not replicated. Nothing was tuned to make it pass.

| Check | How it is measured | Result (main phenotype) |
|---|---|---|
| Has structure | Silhouette vs. a null of 100 column-permuted datasets (must beat the 95th percentile) | 0.36 vs. 0.11 |
| Is stable | Per-cluster bootstrap Jaccard, 200 resamples (≥ 0.75) | 0.92 |
| Travels | Leave-one-cohort-out: medoid transfer vs. native clustering (ARI, Hungarian-matched Jaccard) | ARI 0.50–0.73 |
| Evolves differently | Mixed model `DLCO ~ phenotype × log(time)` | −11.5 points at 12 months |
| Recognisable with few variables | 1–2 question decision tree validated on unseen cohorts | "FVC at 3 m ≤ 83 %?" matches 95–97 % |
| Makes clinical sense | Profiles built from variables **not** used to define the phenotype | Pending review by the clinical team |

### Why this was hard
- **90.5 % of the cells in the merged table are empty.** Each registry used its own case report form, and only **13 clinical variables are shared** by all four cohorts.
- **Who gets a DLCO depends on the hospital:** 36 of 76 CIBERESUCICOVID centres never measure it.
- **Cohorts measure at different times and on different populations**, so we used real dates and fixed clinical ranges instead of per-cohort z-scores (which would erase real differences).
- **Dropout is informative:** patients who don't come back at one year started with a DLCO 6–7 points higher.
- **437 patients appear in two registries**, with the same test recorded twice. Each patient is assigned to a single analysis cohort so that "replication" never reuses the same people.

### Methods at a glance
- **Data cleaning:** 14 documented rules (R01–R14). Every rule logs how many records it touched, and nothing is deleted silently.
- **Clustering:** Gower distance (mixed data, pairwise missingness, no imputation) with fixed clinical ranges and equal weight per domain, then PAM / k-medoids (FasterPAM). The number of clusters `k` is chosen with a permutation null test.
- **Replication:** leave-one-cohort-out, Hungarian matching, ARI and Jaccard with bootstrap confidence intervals.
- **Trajectories:** linear mixed models with random intercept and slope (`statsmodels`), inverse probability weighting for dropout, regression-to-the-mean check and phenotype transitions.
- **Prediction:** Explainable Boosting Machines and L1 logistic regression, validated on unseen cohorts or centres. Reported with ROC-AUC, PR-AUC, calibration (Brier score, slope) and decision curves. Leakage is guarded by code assertions.

---

## Repository structure

> The code, notebooks and full report are written in Spanish.

```
├── main.ipynb                          ← executable end-to-end summary (the notebook presented to the jury)
├── Informe.md                          ← full report: data, methods, results, limitations, jury Q&A
├── config.yaml                         ← every parameter (thresholds, windows, k, seeds, paths)
├── Workflow/
│   ├── 01_analisis_exploratorio.ipynb  ← exploratory analysis: what data, from whom, when, what quality
│   ├── 02_limpieza.ipynb               ← 14 cleaning rules, patient flow (CONSORT), clean tables
│   ├── 03_fenotipado_ciberes.ipynb     ← symptom + function phenotypes in CIBERESUCICOVID (do not travel)
│   ├── 04_fenotipado_comun.ipynb       ← phenotypes using all 3 cohorts together: the main phenotype
│   ├── 05_trayectorias.ipynb           ← DLCO trajectories by phenotype (mixed model, dropout, transitions)
│   ├── 06_modelo_alta.ipynb            ← can the phenotype be predicted at discharge? (EBM, 3 cohorts)
│   ├── 07_modelo_fase_aguda.ipynb      ← adding acute-phase ICU data (EBM, CIBERESUCICOVID only)
│   └── 08_modelo_primera_visita.ipynb  ← at the 3-month visit: who is still impaired at one year?
└── src/                                ← reusable, documented functions
    ├── carga.py                        data loading and long-format follow-up tables
    ├── limpieza.py                     cleaning rules R01–R14 and cleaning log
    ├── fenotipado.py                   Gower distance, PAM, null test, stability
    ├── fenotipado_perfil.py            phenotyping driven by a profile in config.yaml
    ├── pasaporte.py                    cross-cohort replication (Hungarian matching, Jaccard, ARI)
    ├── trayectorias.py                 mixed models, contrasts, IPW
    ├── fichas.py                       phenotype comparison tables (SMD and interpretation)
    └── privacidad.py                   small-cell suppression (N < 10)
```

---

## Principles

- **Reproducible:** every parameter lives in `config.yaml`, with a fixed seed (2026). There are no magic numbers in the code.
- **Replication before discovery:** a phenotype only counts if it passes the passport, and failures are reported too.
- **No imputation of structural missingness:** a variable that a cohort never collects is never invented.
- **Interpretable by design:** real-patient medoids, additive models with one curve per variable, and a one-line clinical rule when that is all you need.
- **Privacy first:** only aggregates are ever shown, with small cells suppressed (`src/privacidad.py`).

## Limitations

- Only 13 variables are comparable across all cohorts, so the replicable phenotype is respiratory only.
- Dropout is corrected under a missing-at-random assumption. Dropout that depends on unobserved future DLCO cannot be ruled out.
- TENACITY is small, so its estimates have wide intervals.
- A fixed 80 % threshold is used instead of the lower limit of normal. Near the threshold, part of the "recovery" may be test variability.
- At the 3-month visit, patients are ranked correctly across cohorts but absolute risk is not: probabilities should be **recalibrated locally** before clinical use.

See [`Informe.md`](Informe.md) (Spanish) for the full report.

## Team

Developed during the Respira Hackathon, organised by **CIBERES** and **AstraZeneca**, by:

- Carlos Palazón Domingo
- Ferran Òdena
- Nil Muriach López


