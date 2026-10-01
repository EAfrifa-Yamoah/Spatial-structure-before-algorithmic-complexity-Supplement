Overview

This repository accompanies our systematic synthesis of evidence on species distribution modelling (SDM) under data-limited marine conditions. The review focuses on the modelling challenges commonly encountered in marine and coastal systems, where observations are often sparse, spatially clustered, and irregularly sampled.

We screened 330 published records and retained 95 studies that met three criteria: (i) a peer-reviewed marine application, (ii) explicit treatment of spatial autocorrelation, and (iii) an available software implementation. Quantitative evidence was extracted from studies reporting direct method comparisons, validation optimism assessments, or sample-size experiments, while the remaining studies informed a structured narrative synthesis.

Key Findings

The review identified four consistent themes:

1. Explicit spatial modelling improves predictive performance more reliably than changing algorithm family.
2. Spatial coverage of sampling is more important than sample size alone, with clustered observations limiting model transferability and predictive reliability.
3. Validation strategies should align with the intended inference goal. Random cross-validation is appropriate for interpolation within sampled conditions but can substantially overestimate performance when models are applied beyond the sampled domain.
4. Minimum data requirements are often modest and method-dependent, with approximately 50–100 well-distributed observations frequently sufficient for tree-based classification approaches.

Purpose of this Repository

This repository contains the supplementary materials supporting the review, including:

1. Data extraction and synthesis files
2. Analysis scripts and reproducible workflows
3. Supporting figures and tables
4. Methodological documentation
5. Additional results not included in the main manuscript
6. Main Conclusion

Across the reviewed literature, predictive reliability under data limitation depends more on how spatial structure and sampling design are handled than on algorithmic sophistication. The repository supports the study's proposed decision framework linking data characteristics, inference objectives, and available sample sizes to appropriate modelling and validation choices.
