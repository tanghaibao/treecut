"""
Draws vertically presented hierarchical tree, along with the associated values
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def _clear_ax(ax):
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    ax.set_axis_off()


class Dendrogram:

    def __init__(self, ext_tree, **kwargs):

        self.tree = ext_tree
        values = self.tree.values

        fig = plt.figure(1, (8,8))
        root = fig.add_axes([0,0,1,1])

        self.draw_tree()
        self.draw_values()
        self.draw_modules()

        _clear_ax(root)


    def _draw_node(self):
        pass


    def _draw_edge(self):
        pass
    

    def draw_tree(self):

        tree_ax = fig.add_axes
        farthest, max_dist = self.tree.node.get_farthest_node()

        # scale the tree so that the height is about .9
        scale = .9 / max_dist

        _clear_ax(tree_ax)

    
    def draw_values(self):

        return None

    
    def draw_modules(self):

        return None


    def savefig(self, image_name, **kwargs):
        plt.savefig(image_name, **kwargs)
