from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import AlignIO
import matplotlib.pyplot as plt
import Logger

logger = Logger.setup_logger()

def build(alignment_path, image_extensions=['pdf']):
    # Load the alignment file
    alignment = AlignIO.read(alignment_path, 'fasta')

    # Calculate the distance matrix
    calculator = DistanceCalculator('identity')
    distance_matrix = calculator.get_distance(alignment)

    # Create a DistanceTreeConstructor object
    constructor = DistanceTreeConstructor(calculator)

    # Construct the phylogenetic tree using the Neighbor-Joining method
    tree = constructor.build_tree(alignment)

    # Save the tree to a file in the specified image formats
    Phylo.draw(tree, do_show=False)
    for ext in image_extensions:
        logger.info(f"Saving phylogenetic tree as 'phylogenetic_tree.{ext}'\n")
        plt.savefig(f'phylogenetic_tree.{ext}')

    # Save the tree in Newick format
    logger.info("Saving phylogenetic tree in Newick format as 'phylogenetic_tree.newick'\n")
    Phylo.write([tree], 'phylogenetic_tree.newick', 'newick')
