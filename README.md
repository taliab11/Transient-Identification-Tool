# Transient-Identification-Tool

## Overview
This repository contains a tool for identifying transients, pulses, and temporary changes over time.
Although originally designed for biological research, specifically for gene-expression analysis, the tool can be adapted for many other fields. 
It identifies genes or other vectorial changes over time that start at a certain level, increase or decrease, and return to a steady level—forming a "mountain" or "valley" shape.

## What is a Transient?
According to its [dictionary definition](https://www.dictionary.com/browse/transient), a "Transient" is something that lasts for only a short time, temporary.
In biology, many processes are transient, including those occurring during cell state transitions, which involve temporary changes in gene expression. Identifying genes with such transient behavior can help us better understand their roles and effects on cellular processes.

## What are the optional algorithms and what are their differences?
The tool offers two algorithms to identify transients, both based on aligning the "transient candidate" (e.g., expression levels of a single gene over time) with a linear line connecting its two ends. The premise is that the less successfully the candidate aligns with the line (i.e., the higher its absolute score), the more likely it is to represent a transient.

### The Euclidean algorithm:
This method measures the sum of all vertical distances between the candidate and the linear line (along the Y-axis, representing gene expression levels) at each time point (X-axis). Essentially, this algorithm detects "non-linearity" and may capture candidates (genes) with upward or downward trends that are not linear but are also not actual transients. Additionally, it may be subject to timeline-scale biases.
### The DTW algorithm:
"Dynamic Time Warping" [(DTW)](https://rtavenar.github.io/blog/dtw.html) is a similarity measure for time series that minimizes the Euclidean distance between aligned series by allowing flexible temporal alignment. Unlike the Euclidean method, DTW identifies "non-linear trends," making it less prone to false positives and potentially more accurate for detecting transients. However, it may also be stricter, missing some "true transients" compared to the Euclidean method.

<p align="center">
  <img src="euc_vs_dtw.png" width="500" style="border-radius: 15px;">
</p>

In the input parameters, you can choose your preferred method (as explained above). If unsure, it is recommended to run the tool twice using both algorithms.

## How is significance measured?
Significance is determined by identifying the "transient candidates" (e.g., genes) that score as significant according to the chosen method.
Using the "Monte Carlo" model, the candidate data is permuted thousands of times, and scores are calculated for each permutation with the chosen algorithm. The score of the actual candidate is then compared to the distribution of scores from these permutations. A candidate is considered significant, with a high probability of being transient, if its P-value is below 0.05 (after FDR adjustment for multiple hypotheses).

## User's input
* The data - A text file formatted as a table (data frame), where rows represent the different "transient candidates" (genes) and columns represent the samples, ordered by time.
* Algorithm choice - Either Euclidean or DTW (as explained above).
* "Monte Carlo" - The desired number of permutations (The default is 5000). A higher number increases accuracy but also runtime.
* Multiple testing correction method (The default is FDR-BH).
* Timestamps - A numeric vector of the ordered time stamps of the samples (in the same time units).
* Repeats columns - The column numbers representing each repeat in the data. The code was built for 2 repeats but can be updated easily for more.
* "Candidate ID" columns - the number of the column in the df containing the "candidate id" (The "gene-symbol" for example).
* Grid name - The desired name for the png file of the grid of plots.

## Output
* The data - The input text file with additional columns for calculated P-values and adjusted P-values.
* plots - A grid of plots, in a png file, visualizing the changes over time and fold-changes for each significant candidate (gene). When there aren't any significant candidates, an appropiate message will be printed and the program will end.

## How to run and required packages
How to run - **ADD

The required packages for running the code are:
* numpy
* pandas
* matplotlib
* seaborn
* dtaidistance
* statsmodels

To install them, prior to running the code, download the requirements.txt file and run:
```
pip install -r requirements.txt
```

## Acknowledgments
* This tool was developed for a project in [Dr. Yaron Antebi's lab](https://www.weizmann.ac.il/molgen/Antebi/).
* This repository was created as part of a [Python programming course](https://github.com/szabgab/wis-python-course-2024-11?tab=readme-ov-file) at WIS.
