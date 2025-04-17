import polars as pl
import matplotlib as mpl
import matplotlib.pyplot as plt

# Read the TSV file into a DataFrame
df = pl.read_csv('data_processed/Bayescode_merged_cov.tsv', separator='\t')
print(f"Number of rows in the Bayecode DataFrame: {df.height}")
# Create a 2x2 grid of subplots
fig, axs = plt.subplots(2, 2, figsize=(12, 10))

# Define columns to plot
columns = [
    "Omega_Aquatic_adaptation_cor",  # Top left
    "Omega_Aquatic_adaptation_partcor",  # Top right
    "Omega_Aquatic_adaptation_ppos_cor",  # Bottom left
    "Omega_Aquatic_adaptation_ppos_partcor"  # Bottom right
]

titles = [
    "Correlation coefficient",
    "Partial correlation coefficient",
    "Posterior probability of positive correlation",
    "Posterior probability of positive partial correlation"
]

# Plot histograms
for i, (column, title) in enumerate(zip(columns, titles)):
    row, col = i // 2, i % 2
    ax = axs[row, col]

    # Determine appropriate x-axis limits and vertical line
    if "ppos" in column:  # Posterior probability
        x_min, x_max = 0, 1
        percentage = df.filter(pl.col(column) > 0.5).height / df.height * 100
        stat_text = f"{percentage:.2f}% of values > 0.5"
    else:  # Correlation coefficient
        x_min, x_max = -1, 1
        vertical_line = 0
        percentage = df.filter(pl.col(column) > 0).height / df.height * 100
        stat_text = f"{percentage:.2f}% of values > 0"
        ax.axvline(vertical_line, color='black', linestyle='-', linewidth=2)

    # Plot histogram
    ax.hist(df[column].to_numpy(), bins=50, edgecolor='black')

    # Set limits and labels
    ax.set_xlim(x_min, x_max)
    ax.set_xlabel(title)
    ax.set_ylabel('Frequency')
    ax.set_title(f"{stat_text}")

# Add column labels
fig.text(0.30, 0.95, 'Correlation', ha='center', fontsize=14, fontweight='bold')
fig.text(0.775, 0.95, 'Partial Correlation', ha='center', fontsize=14, fontweight='bold')

# Add row labels
fig.text(0.02, 0.725, 'Correlation coefficients', va='center', rotation=90, fontsize=14, fontweight='bold')
fig.text(0.02, 0.275, 'Posterior probabilities', va='center', rotation=90, fontsize=14, fontweight='bold')

plt.tight_layout(rect=[0.05, 0.05, 1, 0.95])  # Adjust layout to make room for labels
plt.savefig('results/Bayescode_results.pdf', dpi=300, bbox_inches='tight')
plt.show()

df_hyphy = pl.read_csv('data_processed/Hyphy_merged.tsv', separator='\t')
print(f"Number of rows in the Hyphy DataFrame: {df_hyphy.height}")
# Assert that the number of rows in hyphy DataFrame is the same as the number of rows in the BayesCode DataFrame
assert df_hyphy.height == df.height, "The number of rows in the Hyphy DataFrame does not match the BayesCode DataFrame."

# Create a 1x2 grid of subplots
fig, axs = plt.subplots(2, 2, figsize=(12, 10))
# Histogram of p-value distributions and k
ax = axs[0, 0]
ax.hist(df_hyphy["p-value"].to_numpy(), bins=50, edgecolor='black')
ax.set_xlim(0, 1.0)
ax.set_xlabel("p-value from Relax")
ax.set_ylabel('Frequency')

ax = axs[0, 1]
ax.hist(df_hyphy["k"].to_numpy(), bins=50, edgecolor='black')
ax.set_xlabel("k parameter from Relax")
ax.set_ylabel('Frequency')


# Merge the two DataFrames on the "id" column
# Plot posterior probability from OrthoMam as function of p-value from hyphy as a density plot
relax_column = "p-value"
relax_title = "p-value from Relax"
for i, (title, column) in enumerate([("Correlation", "Omega_Aquatic_adaptation_ppos_cor"),
                                     ("Partial correlation", "Omega_Aquatic_adaptation_ppos_partcor")]):
    ax = axs[1, i]
    # Use a density plot to visualize the density of points, use a logarithmic scale for the color map
    ax.set_xlim(0, 1.0)
    ax.set_ylim(0.0, 1.0)
    merged_df = df_hyphy.join(df, on="id", how="inner")
    x = merged_df[relax_column].to_numpy()
    y = merged_df[column].to_numpy()
    ax.hist2d(x, y, bins=10, norm=mpl.colors.LogNorm(), cmap="Blues")
    ax.set_xlabel(relax_title)
    ax.set_ylabel(f"Posterior probability from BayesCode")
    ax.set_title(title)

plt.savefig('results/Hyphy_results.pdf', dpi=300, bbox_inches='tight')
plt.show()
