# ACreceptor

**ACreceptor** is an R package designed to assist researchers in inferring receptor activation and cellular communication. The package provides efficient tools for data processing, correlation analysis, and visualization, enabling users to explore the relationships between receptor activity and cell-cell communication.

## Features
- Infers receptor activation and its role in cell-cell communication.
- Supports efficient processing and analysis of receptor-related data.
- Provides visualization tools to help users interpret receptor activity and communication patterns.

## Installation

### Install from GitHub using devtools
To install the development version of **ACreceptor** from GitHub, you can use the `devtools` package:
```r
# Install devtools package (if not already installed)
install.packages("devtools")

# Install ACreceptor from GitHub
devtools::install_github("xuansuyang/ACreceptor")

## Usage
After installing the package, you can load it and start using the available functions. Here’s a simple example:
library(ACreceptor)
model=readRDS(system.file("extdata", "human_model.rds", package = "ACreceptor"))
expr=read.table("input/counts/GSE180698_counts.txt",header = T,row.names = 1)
expr=cpm(expr)
expr=cpm_normalized(expr)
ac=activity(expr,model)
