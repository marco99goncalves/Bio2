from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo
import matplotlib.pyplot as plt
import os
import Logger

logger = Logger.setup_logger()

def make_pretty_tree(tree):
    """
    Make adjustments to the tree to improve its aesthetics.
    
    :param tree: The Phylo tree object to be prettified.
    """
    # Set branch line colors and widths
    for clade in tree.find_clades():
        # Replace underscores with spaces, capitalize the first word, and make the rest lowercase
        words = clade.name.replace('_', ' ').split()
        words[0] = words[0].capitalize()
        for i in range(1, len(words)):
            words[i] = words[i].lower()
        clade.name = ' '.join(words)

        if not clade.is_terminal():
            # Internal branches
            clade.color = 'blue'
        else:
            # Terminal branches
            clade.color = 'green'
            
    # Set the style of the tree
    tree.ladderize()  # Ladderize to organize the tree

def build(alignment_path, export_formats=["pdf"], matrix_name="BLOSUM62"):
    # Read the alignment file
    alignment = AlignIO.read(alignment_path, "fasta")
    
    matrix_name = matrix_name.lower()

    # Calculate the distance matrix using the 'blosum62' model
    calculator = DistanceCalculator(matrix_name)
    distance_matrix = calculator.get_distance(alignment)
    
    # Construct the phylogenetic tree using the Neighbor-Joining algorithm
    constructor = DistanceTreeConstructor(calculator, 'nj')
    tree = constructor.build_tree(alignment)
    
    # Make the tree pretty
    make_pretty_tree(tree)
    
    # Set up the Matplotlib figure and axis
    plt.figure(figsize=(10, 10))
    ax = plt.subplot(1, 1, 1)
    
    # Draw the tree with the prettified style settings
    Phylo.draw(tree, do_show=False, axes=ax)
    
    # Set the title and font size for the labels
    ax.set_title('Phylogenetic Tree for the P68871 gene', fontsize=16)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontsize(12)
    
    # Export the tree to the specified file formats
    file_basename = "Outputs/PhylogeneticTree"
    for export_format in export_formats:
        fileName = f"{file_basename}.{export_format}"
        plt.savefig(f"{fileName}")
        logger.info(f"Exported tree to {fileName}")
    
    plt.close()
