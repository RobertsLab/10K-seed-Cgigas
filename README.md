# 10K Seed Hardening Project 

## Overview

This repository contains research data and analysis on Pacific oyster (*Crassostrea gigas*) hardening treatments designed to improve survival and stress tolerance in aquaculture operations.

## Industry Recommendations

**🔗 [VIEW INDUSTRY RECOMMENDATIONS](INDUSTRY_RECOMMENDATIONS.md)** 

**📋 [QUICK REFERENCE PROTOCOL](QUICK_REFERENCE.md)**

For practical guidance on implementing oyster hardening procedures based on this research, see the comprehensive industry recommendations document.

**Quick Summary:** Fresh water + temperature combined hardening shows the most promise for improving heat stress survival, though benefits are modest and should be evaluated against implementation costs.

## Research Components

### Laboratory Studies
- **Survival Analysis** (`scripts/survival.md`) - Kaplan-Meier survival curves and Cox proportional hazards analysis
- **Metabolic Rate Analysis** (`scripts/resazurin.md`) - Resazurin assay to measure metabolic responses to stress
- **Temperature Logger Data** (`scripts/loggers.Rmd`) - Environmental monitoring during experiments

### Field Studies  
- **Outplant Survival** (`scripts/outplant-survival.Rmd`) - Long-term field performance of hardened oysters

### Treatments Tested
- **Control** - No hardening treatment
- **Fresh Water** - Osmotic hardening with fresh water exposure
- **Temperature** - Thermal hardening with elevated temperature exposure  
- **Fresh Water + Temperature** - Combined osmotic and thermal hardening
- **Immune** - Immune system stimulation

## Key Findings

- No significant survival differences between treatments under normal conditions (18°C)
- Significant differences under heat stress (42°C), with fresh water + temperature showing best performance (15.5% vs 5.5% for temperature-only)
- Metabolic rates can predict mortality risk with high accuracy (AUC = 0.996)
- Long-term field benefits are limited, suggesting hardening effects are most relevant for acute stress events

## Data Structure

```
data/
├── resazurin/          # Metabolic rate assay data
├── survival/           # Laboratory survival experiment data  
├── outplant/           # Field deployment survival data
├── loggers/            # Environmental monitoring data
└── sample-metadata/    # Sample tracking and treatment assignments
```

## Usage

For researchers interested in the detailed methodology and statistical analysis, examine the R Markdown files in the `scripts/` directory. For industry applications, start with the [Industry Recommendations](INDUSTRY_RECOMMENDATIONS.md) document.
