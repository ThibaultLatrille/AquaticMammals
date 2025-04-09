import polars as pl
import matplotlib.pyplot as plt

# Read the TSV file into a DataFrame
df = pl.read_csv('data_processed/OrthoMam_merged_cov.tsv', separator='\t')

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
plt.savefig('results/correlation_histograms.pdf', dpi=300, bbox_inches='tight')

# Save the name of the CDS with Omega_Aquatic_adaptation_ppos_partcor > p
p = 0.999  # Set the threshold for filtering
filtered_df = df.filter(pl.col("Omega_Aquatic_adaptation_ppos_partcor") > p)
print(f"Number of CDS with Omega_Aquatic_adaptation_ppos_partcor > {p}: {filtered_df.height}")
filtered_df.select("id").write_csv("results/filtered_CDS.csv")
plt.show()