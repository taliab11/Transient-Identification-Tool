#Loading required libraries:
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from dtaidistance import dtw
from statsmodels.stats.multitest import multipletests
import argparse

#"Eucledean" method:
#Creating "GeneScore" function - receives a vector of the candidate levels over time and retrieves a score - a sum of the distances from the linear line between the two ends of the distribution. values above the line will have positive distances, and below negative distances.
def cand_score(cand_over_time):
    cand_over_time = np.array(cand_over_time)
    mean_start = np.mean(cand_over_time[:len(cand_over_time) // 5])
    mean_end = np.mean(cand_over_time[-(len(cand_over_time) // 6):])
    slope = (mean_end - mean_start) / len(cand_over_time)
    intercept = mean_start - slope
    x = np.arange(1, len(cand_over_time) + 1)
    distances = cand_over_time - (slope * x + intercept)
    return np.sum(distances)

#Creating "LinearReference" function for the DTW method and the FC calculations - receives a vector of the candidate levels over time and retrieves the corresponding y values of a linear line between the two ends of the distribution:
def linear_reference(cand_over_time):
    cand_over_time = np.array(cand_over_time)
    mean_start = np.mean(cand_over_time[:len(cand_over_time) // 5])
    mean_end = np.mean(cand_over_time[-(len(cand_over_time) // 6):])
    slope = (mean_end - mean_start) / len(cand_over_time)
    intercept = mean_start - slope
    x = np.arange(1, len(cand_over_time) + 1)
    return slope * x + intercept

#"DTW" method:
#Creating "DTWscore" function - receives a vector of the candidate levels over time and retrieves the min distance between the vector and a linear line between the two end of the distribution, using the DTW algorithm:
def dtw_score(cand_over_time):
    query = np.array(cand_over_time)
    reference = linear_reference(query)
    return dtw.distance(query, reference)

#"Monte Carlo" (Permutations) - calculating mixing options scores - used for both methods:
def mix_scores(cand_over_time, method, monte_carlo=5000):
    scores = []
    for _ in range(monte_carlo):
        shuffled = np.random.permutation(cand_over_time)
        if method == "Euclidean":
            scores.append(cand_score(shuffled))
        elif method == "DTW":
            scores.append(dtw_score(shuffled))
    return np.array(scores)

#Calculating the p-value (two-sided) - used for both methods:
def score_pvalue(cand_over_time, method, monte_carlo=5000):
    observed_score = cand_score(cand_over_time) if method == "Euclidean" else dtw_score(cand_over_time)
    null_distribution = mix_scores(cand_over_time, method, monte_carlo)
    return np.mean(np.abs(null_distribution) >= np.abs(observed_score))

# The "main" analysis function:
# A general function that sends the data (in a df, ordered by time) to the chosen calculation measure and returns the df with a Pvalue and adjusted Pvalue (according to the chosen method) columns:
def transient_pvalue(df, method, monte_carlo=5000, adj_method='fdr_bh', repeat1_cols=[], repeat2_cols=[]):
    pvalues = df.apply(lambda row: score_pvalue(row[repeat1_cols + repeat2_cols], method, monte_carlo), axis=1)
    adj_pvalues = multipletests(pvalues, method=adj_method)[1]
    df["Pvalue"] = pvalues
    df["adj_Pvalue"] = adj_pvalues
    return df

#calculating "Fold Change" for the plots by the mean of the two highest FC points (by proportion) between the candidate levels and the linear line connecting the two ends:
def max_fc(cand_over_time):
    cand_over_time = np.array(cand_over_time) + 1  # Shift to avoid log issues
    ref_line = linear_reference(cand_over_time)
    distances = cand_over_time - ref_line
    max_indices = np.argsort(np.abs(distances))[-2:]
    fc_values = cand_over_time[max_indices] / ref_line[max_indices]
    return np.mean(fc_values)

# Visualization of the results:
#creating a grid of plots of the significant P-Value genes in decreasing log2(FC) order - including the results df, the x-axis time values (as integers in the same unit), the columns for each repeat and the "candidate id" column. The grid will be saved as a png file in the code's folder:
#(The function was built for 2 repeats but can be updated easily for more)
def plot_grid(df, x_axis, repeat1_cols, repeat2_cols, candidate_id_col, grid_name):
    df = df.sort_values(by='adj_Pvalue', ascending=True)

    # to choose the dimensions of the grid we will find the "whole square root" that is equal or higher than the sqrt of the sum of the significant Pvalues:
    significant_cands = df[df["adj_Pvalue"] <= 0.05].head(int(np.ceil(np.sqrt(len(df)))))
    
    # Sort genes by fold change in decreasing order
    sorted_cands = significant_cands.index[np.argsort(-np.log2(significant_cands.apply(lambda row: max_fc(row[repeat1_cols + repeat2_cols]), axis=1)))]
    if len(sorted_cands) == 0:
        print("No significant transient candidates found for plotting.")
        return

    n = int(np.ceil(np.sqrt(len(sorted_cands))))
    fig, axes = plt.subplots(n, n, figsize=(15, 15))
    if isinstance(axes, np.ndarray):
        axes = axes.flatten()
    else:
        axes = [axes]
    
    for i, cand in enumerate(sorted_cands):
        ax = axes[i]
        time_points = np.log2(x_axis) #log on X-axis (for visualization)
        y_repeat1 = df.iloc[df.index[cand], repeat1_cols].values
        y_repeat2 = df.iloc[df.index[cand], repeat2_cols].values
        fc_value = max_fc(df.iloc[df.index[cand], repeat1_cols + repeat2_cols])
        
        ax.scatter(time_points, y_repeat1, color='green', label='Repeat1') # Define colors for repeats
        ax.scatter(time_points, y_repeat2, color='orange', label='Repeat2')
        ax.plot(time_points, y_repeat1, color='green', alpha=0.7)
        ax.plot(time_points, y_repeat2, color='orange', alpha=0.7)
        ax.set_title(f"{df.iloc[df.index[cand], candidate_id_col]} (log2(FC): {fc_value:.3f})")
        ax.set_xticks([])
        ax.set_yticks([])
    
    plt.tight_layout()
    plt.savefig(grid_name)
    plt.close()

# The function that runs the analysis and visualization:

# The parameters: df = the data, ordered by time, method = the chosen calculation measure ("Euclidean"/"DTW"), MonteCarlo = the number of permutation of each gene for Pvalue calculation, adj_method = chosen adj. method for "multiple hypothesis" (can receive all the "multipletests" optional methods), x_axis = the x-axis time values (as integers in the same unit), repeat1_cols = the column of rep1 in the df, repeat2_cols = the column of rep2 in the df, candidate_id_col = the number of the column in the df containing the "candidate id".
def transient_analysis(df, method, monte_carlo=5000, adj_method='fdr_bh', x_axis=[], repeat1_cols=[], repeat2_cols=[], candidate_id_col=0, grid_name='plot_grid.png'):
    updated_df = transient_pvalue(df, method, monte_carlo, adj_method, repeat1_cols, repeat2_cols)
    print(updated_df)
    plot_grid(updated_df, x_axis, repeat1_cols, repeat2_cols, candidate_id_col, grid_name)
    return updated_df

parser = argparse.ArgumentParser()
parser.add_argument('--df', required=True, type=str)
parser.add_argument('--algorithm', required=True, type=str)
parser.add_argument('--monte_carlo', required=False, type=int, default=5000)
parser.add_argument('--adj_method', required=False, type=str, default="fdr_bh")
parser.add_argument('--time_stamps', required=True, type=int, nargs="+")
parser.add_argument('--repeat1_cols', required=True, type=int, nargs="+")
parser.add_argument('--repeat2_cols', required=False, type=int, nargs="+", default=[])
parser.add_argument('--candidate_id_col', required=False, type=int, default=0)
parser.add_argument('--grid_name', required=False, type=str, default="plot_grid")

args = parser.parse_args()
df = pd.read_csv(args.df, sep="\t", header=0)
method = args.algorithm
monte_carlo = args.monte_carlo
adj_method = args.adj_method
x_axis = args.time_stamps
repeat1_cols = args.repeat1_cols
repeat2_cols = args.repeat2_cols
candidate_id_col = args.candidate_id_col
grid_name = args.grid_name + ".png"

data_pvalue = transient_analysis(df, method, monte_carlo, adj_method, x_axis, repeat1_cols, repeat2_cols, candidate_id_col, grid_name)