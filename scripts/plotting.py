import logging
import os
from typing import Dict

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from colours import superkingdom_colors, phylum_colors, crassvirales_color
from logging_utils import setup_logging
from utils import time_it

logger = logging.getLogger()  # Get the configured logger


def plot_bacterial_ratios_vs_threshold(concatenated_table: str, output_dir: str, tree_type: str) -> None:
    """Plot bacterial ratios vs thresholds and save the figure."""
    
    figures_dir = os.path.join(output_dir, 'figures')
    os.makedirs(figures_dir, exist_ok=True)

    logger.info("Loading concatenated table from %s", concatenated_table)
    
    try:
        df = pd.read_csv(concatenated_table, sep='\t', index_col=False)
        logger.info(f"First rows of df:\n{df.head()}\n")
    except Exception as e:
        logger.error("Error loading concatenated table: %s", str(e), exc_info=True)
        return
    
    colors = {
        'ratio_Bacteroidetes_to_total': phylum_colors['Bacteroidetes'],  
        'ratio_Actinobacteria_to_total': phylum_colors['Actinobacteria'],  
        'ratio_Bacillota_to_total': phylum_colors['Bacillota'],  
        'ratio_Proteobacteria_to_total': phylum_colors['Proteobacteria'],  
        'ratio_Other_to_total': 'gray',  
        'ratio_bacterial_to_total': superkingdom_colors['Bacteria'],  
        'crassvirales_ratio': crassvirales_color,  
        'ratio_viral_to_total': superkingdom_colors['Viruses']  
    }

    df_long = pd.melt(df, id_vars=['threshold'], value_vars=list(colors.keys()),
                      var_name='Type', value_name='Ratio (%)')
    
    logger.debug("First rows of df_long before resetting index:\n%s", df_long.head())
    logger.debug("Data types:\n%s", df_long.dtypes)
    
    df_long = df_long.reset_index(drop=True)
    
    
    # logger.info("First rows of df_long:\n%s", df_long.head().to_string())
    logger.debug("First rows of df_long after resetting index:\n%s", df_long.head())
    logger.debug("Data types:\n%s", df_long.dtypes)

    plt.figure(figsize=(12, 8))
    
    try:
        logger.debug(f"df_long['threshold'] shape: {df_long['threshold'].shape}, ndim: {df_long['threshold'].ndim}")
        logger.debug(f"df_long['Ratio (%)'] shape: {df_long['Ratio (%)'].shape}, ndim: {df_long['Ratio (%)'].ndim}")

        # sns.lineplot(x=df_long["threshold"].values.flatten(), y=df_long["Ratio (%)"].values.flatten(), hue='Type', data=df_long, palette=colors)
        sns.lineplot(x=df_long["threshold"], y=df_long["Ratio (%)"], hue='Type', data=df_long, palette=colors)

    except Exception as e:
        logger.error("Error while generating plot: %s", str(e), exc_info=True)
        return

    plt.title(f'Bacterial, Viral, and Crassvirales Ratios vs Thresholds ({tree_type.capitalize()} Tree)')
    plt.xlabel('Threshold of Crassvirales proteins in the clade (%)')
    plt.ylabel('Ratio to Total Members (%)')
    plt.legend(title='Type', loc='best')

    output_file = os.path.join(figures_dir, f'bacterial_viral_crassvirales_ratios_vs_threshold_{tree_type}.png')
    
    try:
        plt.savefig(output_file, dpi=300)
        plt.close()
        logger.info("Plot saved to %s", output_file)
    except Exception as e:
        logger.error("Error saving plot: %s", str(e), exc_info=True)


def plot_crassvirales_bacterial_viral_ratios_vs_threshold(concatenated_table: str, output_dir: str,
                                                          tree_type: str) -> None:
    """Plot Crassvirales, bacterial, and viral ratios vs thresholds and save the figure."""
    # Create the figures directory if it doesn't exist
    figures_dir = os.path.join(output_dir, 'figures')
    os.makedirs(figures_dir, exist_ok=True)

    # Load the concatenated table
    df = pd.read_csv(concatenated_table, sep='\t')

    # Define the dictionary of colors for the lines
    colors: Dict[str, str] = {
        'crassvirales_ratio': crassvirales_color,  # Blue
        'ratio_bacterial_to_total': superkingdom_colors['Bacteria'],  # Orange
        'ratio_viral_to_total': superkingdom_colors['Viruses'],  # Green
        'ratio_other_to_total': 'gray'
    }

    # Prepare the DataFrame in long format for seaborn plotting
    df_long = pd.melt(df, id_vars=['threshold'], value_vars=list(colors.keys()),
                      var_name='Ratio Type', value_name='Ratio (%)')

    # Plot using seaborn
    plt.figure(figsize=(10, 6))
    sns.lineplot(x='threshold', y='Ratio (%)', hue='Ratio Type',
                 data=df_long, palette=colors)

    # Customize the plot
    plt.title(f'Crassvirales, Bacterial, Viral, and Other Ratios vs Thresholds ({tree_type.capitalize()} Tree)')
    plt.xlabel('Threshold of Crassvirales proteins in the clade (%)')
    plt.ylabel('Ratio to Total Members (%)')
    plt.legend(title='Ratio Type', loc='best')

    # Save the figure
    output_file = os.path.join(figures_dir, f'crassvirales_bacterial_viral_ratios_vs_threshold_{tree_type}.png')
    plt.savefig(output_file, dpi=300)
    plt.close()

    # print(f"Plot saved to {output_file}")


def plot_number_of_clades_vs_threshold(concatenated_table: str, output_dir: str, tree_type: str) -> None:
    """Plot the number of clades found vs thresholds and save the figure."""
    figures_dir = os.path.join(output_dir, 'figures')
    os.makedirs(figures_dir, exist_ok=True)

    df = pd.read_csv(concatenated_table, sep='\t')

    # Count the number of clades for each threshold
    clade_counts = df.groupby('threshold').size().reset_index(name='Number of Clades')

    plt.figure(figsize=(10, 6))
    sns.lineplot(x='threshold', y='Number of Clades', data=clade_counts, marker='o')

    plt.title(f'Number of Clades Found vs Thresholds ({tree_type.capitalize()} Tree)')
    plt.xlabel('Threshold of Crassvirales proteins in the clade (%)')
    plt.ylabel('Number of Clades')

    output_file = os.path.join(figures_dir, f'number_of_clades_vs_threshold_{tree_type}.png')
    plt.savefig(output_file, dpi=300)
    plt.close()

    # print(f"Plot saved to {output_file}")


def plot_number_of_members_boxplot(concatenated_table: str, output_dir: str, tree_type: str) -> None:
    """Plot boxplots of number_of_members vs thresholds and save the figure."""
    figures_dir = os.path.join(output_dir, 'figures')
    os.makedirs(figures_dir, exist_ok=True)

    # Load the concatenated table
    df = pd.read_csv(concatenated_table, sep='\t')

    plt.figure(figsize=(12, 8))
    sns.boxplot(x='threshold', y='number_of_members', data=df)

    plt.title(f'Boxplot of Number of Members vs Thresholds ({tree_type.capitalize()} Tree)')
    plt.xlabel('Threshold of Crassvirales proteins in the clade (%)')
    plt.ylabel('Number of Members in the clade')

    output_file = os.path.join(figures_dir, f'number_of_members_boxplot_{tree_type}.png')
    plt.savefig(output_file, dpi=300)
    plt.close()


@time_it("Generating plots")
def generate_plots(output_paths: Dict[str, str], tree_type: str) -> None:
    """Generate and save all relevant plots.

    Args:
        output_paths (Dict[str, str]): Dictionary containing output paths for the various files.
        tree_type (str): The type of the tree being analyzed (e.g., 'rooted', 'unrooted', 'midpoint').
    """
    concatenated_table = output_paths['biggest_non_intersecting_clades_all']

    plot_bacterial_ratios_vs_threshold(concatenated_table, output_paths['output_dir'], tree_type)
    plot_crassvirales_bacterial_viral_ratios_vs_threshold(concatenated_table, output_paths['output_dir'], tree_type)
    plot_number_of_clades_vs_threshold(concatenated_table, output_paths['output_dir'], tree_type)
    plot_number_of_members_boxplot(concatenated_table, output_paths['output_dir'], tree_type)
