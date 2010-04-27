"""
Draws vertically presented hierarchical tree, along with the associated values
"""

import itertools

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def clear_ax(ax):
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    ax.set_axis_off()


class Dendrogram(object):

    def __init__(self, ext_tree, **kwargs):

        self.tree = ext_tree
        values = self.tree.values

        fig = plt.figure(1, (8,8))
        root = fig.add_axes([0,0,1,1])
        tree_ax = fig.add_axes([0,.5,1,.5])

        self.draw_tree(tree_ax)
        self.draw_values()
        self.draw_modules()

        clear_ax(root)


    def _draw_node(self):
        pass


    def _draw_edge(self):
        pass
    

    def draw_tree(self, ax):

        t = self.tree.node
        farthest, max_dist = t.get_farthest_node()
        margin = .05
        xstart = margin
        ystart = 1 - margin
        canvas = 1 - 2 * margin
        # scale the tree so that the height is about .9
        scale = canvas / max_dist

        num_leaves = len(t.get_leaf_names())
        xinterval = canvas / (num_leaves - 1)

        i = itertools.count()
        coords = {}
        for n in t.traverse("postorder"):
            dist = n.get_distance(t)
            yy = ystart - scale * dist

            if n.is_leaf():
                xx = xstart + i.next() * xinterval
            else:
                children = [coords[x] for x in n.get_children()]
                children_x, children_y = zip(*children)
                min_x, max_x = min(children_x), max(children_x)
                # plot the horizontal bar
                ax.plot((min_x, max_x), (yy, yy), "g-")
                # plot the vertical bar
                for cx, cy in children:
                    ax.plot((cx, cx), (cy, yy), "r-")
                xx = np.mean(children_x)

            coords[n] = (xx, yy)

        clear_ax(ax)

    
    def draw_values(self):

        return None

    
    def draw_modules(self):

        return None


    def savefig(self, image_name, **kwargs):

        plt.savefig(image_name, **kwargs)
