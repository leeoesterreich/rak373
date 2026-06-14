# Selected Representative Longitudinal cfDNA CNV Heatmap

## Purpose

This workflow generates a gene-level copy-number alteration (CNA/CNV) heatmap for selected representative cfDNA cases using longitudinal CNVKit `.annotated.cnr` files.

The goal is to visualize CNV states across serial plasma samples for selected `Patient_Identifier` values while keeping all longitudinal samples from the same patient together.

The heatmap includes:

* gene-level CNV state calls for selected breast cancer genes,
* longitudinal sample ordering by patient and blood draw,
* patient-level separation using dotted vertical lines,
* top annotations for histology, molecular subtype, and tumor fraction (`TFx`).

---

## Representative longitudinal ichorCNA tumor fraction trajectory categories

The selected representative patients were grouped using the following operational rules based on longitudinal ichorCNA tumor fraction (`TFx`) patterns.

| Category                     | Operational rule used                                                                           |
| ---------------------------- | ----------------------------------------------------------------------------------------------- |
| **Strong increase**          | Final TF - baseline TF = 15 percentage points; max TF = 10%                                     |
| **Strong decrease**          | Baseline TF - final TF = 15 percentage points; baseline TF = 10%                                |
| **Transient peak then drop** | Highest TF occurs at an internal blood draw and drops =15 percentage points by final draw       |
| **Rebound / fluctuating**    | At least one direction change or internal trough; TF range =15 percentage points; max TF =10%   |
| **Persistent high**          | Median TF =30%, minimum TF =20%, and minimum-to-median TF ratio =0.70 across serial blood draws |
| **Stable low**               | Maximum TF <10% and TF range <5%                                                                |

---

## Notes

This workflow is designed for visualization and case selection. It should not be used alone to make biological conclusions without checking the original CNVKit output, tumor fraction, treatment history, and clinical timeline.
