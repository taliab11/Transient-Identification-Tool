# Transient-Identification-Tool

## Overview
This repository contains a tool for identification of transients, pulse and temporary changes, over time.  
Although this tool was originally designed for Biology research and specifically for gene-expression analysis, it can be relevant and implemented in many other fields. 
The tool identifies genes, or other vectorial change over time, that begins in a certain level, increases or decreases and returns to a steady level (expressed in a shape of a "mountain" or a "valley").

## What is a Transient?
According to its [dictionary definition](https://www.dictionary.com/browse/transient), a "Transient" is something that lasts for only a short time, temporary.
In Biology, there are many transient processes, including those that accure during cell state transition and aqcuire transient changes in gene expression. Identification of genes that change transiently during these different processes, can help us better understand their role and effect on the cell.

## What are the optional algorithms and what are their differences?
The tool, and both of the optional algorithms within it, are based on alignment of the "transient candidate" (for example, expression levels of a single gene over time) with a linear line connecting its two ends. The rationale is that the less the linear alignment is "successful", its score will be higher in absolute value and the probability of the candidate of being transient arises.
### The Euclidean algorithm:
This method measures the sum of all distances between the candidate and the linear (in parallel to the Y axis, which represents the gene reads in gene-expression analysis) by their time points (the X axis). Basically, this algorithm identifies "non-linearity" and therefore may capture candidates (genes) with an upward or downward trend that aren't linear but also aren't actual transients. Moreover, this algorithm is subject to timeline scale biases.
### The DTW algorithm:
"Dynamic Time Warping" [(DTW)](https://rtavenar.github.io/blog/dtw.html) is a similarity measure between time series which seeks for the temporal alignment that minimizes Euclidean distance between aligned series. Using this algorithm, the exact timestamps at which observations occur are disregarded, only their ordering matters. Unlike the euclidean method, DTW identifies "non-lineary trends" and therefore is less prone to "false-positives" and might be more accurate for identifying transients. It also might be more "strict" and capture less "true transients" compared to the euclidean method.

<p align="center">
  <img src="euc_vs_dtw.png" width="500" style="border-radius: 15px;">
</p>

In the input parameters you'll be able to chose your preferred method (as explained below), when not sure the recommendation is to run it twice, using both algorithms.

## How is the significance measured?
The transients will be the significant candidates (genes) according to the chosen method.
The significance is measured using the "Monte Carlo" model by permutating the candidate (gene) data thousands of times and calculating the score for each permutation, using the same chosen algorithm. The score of the "real" candidate will be compared to the total permutations score distribution. It will be considered significant, with high probability of being transient, with a P-value below 0.05 (after FDR adjustment for multiple hypothsis).

## User's input
* The data - A text file organized as a table (data frame) where the rows are the different "transient candidates" (genes) and the columns are the samples, ordered by time.
* Time-stamps - A vector of the ordered time stamps of the samples.
* The chosen method - Euclidean or DTW (as explained above).
* "Monte Carlo" - The requested number of permutation - the higher the number the longer the running time and the better the accuracy.

## The output
* The data - The input text file with added columns of the calculated P-values and adjusted P-values.
* A grid of plots presenting the change over time and fold-change for each significant candidate (gene).

## How to run and required packages
**ADD

## Acknowledgments
* This tool was created for a project in [Dr Yaron Antebi's lab](https://www.weizmann.ac.il/molgen/Antebi/).
* This repository was created as part of a [Python programming course](https://github.com/szabgab/wis-python-course-2024-11?tab=readme-ov-file) at WIS.
