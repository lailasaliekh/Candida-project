# Standard modules
import matplotlib.pyplot as plt
import numpy as np
from glob import glob
import pandas as pd
import os
import sys

 #Third party
from shapely.geometry import Polygon

# Custom modules
from RodShapedBacteria import RodShapedBacterium
from findCriticalLengthProb import getCells
import utilities as ut
import fastCellPlotting as fcp
from DistributionFunctions import computeColonyContour,getColonyDensity
from generalPlotting import addNematicDirector

ut.setMPL()
#########################ss##############this is plot chained and non chained cells colour coded so no longer by orientation
save_dir="/Users/s2507701/OneDrive - University of Edinburgh/Confined_channels_2025/snapshots/w400/"
# List of specific cell IDs to be colored yellow
yellow_cells_ids = [ 2878,2283 , 3214] 
# def chaining_status_to_color(cell):
#     """
#     Returns red if the cell is not chained (no upper or lower link),
#     and green if it is chained (has either upper or lower link).
#     """
#     upper_link = cell.upper_link
#     lower_link = cell.lower_link
#     if (upper_link is None or np.isnan(upper_link) or upper_link == "None") and (lower_link is None or np.isnan(lower_link) or lower_link == "None"):
#         return (0, 1, 0, 1)  # Red
#     else:
#         return (1, 0, 0, 1)  # Green
    

def plotCells(file, add_director=True):
    RodShapedBacterium.sus_vis_radius_factor = 0.7
    streamplot = (True & add_director)

    fig, ax = plt.subplots(1, 1, figsize=ut.getFigsize(cols=2),
                           constrained_layout=True,
                           facecolor='k')
    
    # Load and inspect data
    dat = pd.read_csv(file, sep='\t')
    print("Data loaded:")
    print(dat.head())

    cells = getCells(file)

    # Confinement bounds (adjust as needed)
    y_min_wall = -0.09
    y_max_wall = 20.1

    # Filter cells to stay within y boundaries
    # cells = [cell for cell in cells if y_min_wall <= cell.pos_y <= y_max_wall]
    print("Cells loaded:", len(cells))

    # for cell in cells:
    #     cell.colour = chaining_status_to_color(cell)
    #     # Debugging: Print the chaining status and assigned color
    #     print(f"Cell ID: {cell.cell_id}, Upper Link: {cell.upper_link}, Lower Link: {cell.lower_link}, Color: {cell.colour}")

    maxx = dat['pos_x'].max() + 3
    minx = dat['pos_x'].min() - 3
    maxy = dat['pos_y'].max() + 3
    miny = dat['pos_y'].min() - 3
    X = maxx - minx
    X_c = 0.5 * (maxx + minx)
    Y = maxy - miny
    Y_c = 0.5 * (maxy + miny)
    if (X >= Y):
        Y = X
        miny = Y_c - 0.5 * Y
        maxy = Y_c + 0.5 * Y
    else:
        X = Y
        minx = X_c - 0.5 * X
        maxx = X_c + 0.5 * X
    fcp.addAllCellsToPlot(cells, ax, ax_rng=maxx - minx, show_id=True, ec='w')
    fcp.addChainLinks(cells, ax, ax_rng=maxx - minx, show_id=True, colour='w')

    # Add markers on top of cells with IDs in yellow_cells_ids
    # for cell in cells:
    #     if cell.cell_id in yellow_cells_ids:
    #         ax.plot(cell.pos_x, cell.pos_y, marker='o', markersize=6, color='yellow', markeredgecolor='black', markeredgewidth=0.2,alpha=0.6)
    scale=1
# Add horizontal walls, extending from minx to maxx
    wall_color = 'white'
    wall_linewidth = 1

    ax.plot([minx, maxx], [-0.09, -0.09], color=wall_color, linewidth=wall_linewidth, alpha=0.6,linestyle='--')
    # ax.plot([minx, maxx], [80.1, 80.1], color=wall_color, linewidth=wall_linewidth, alpha=0.6,linestyle='--')
    # Vertical walls
    # ax.plot([-75.1/scale, -75.1/scale], [-.10/scale, 150/scale], color=wall_color, linewidth=wall_linewidth, alpha=0.6)
    # ax.plot([75.1/scale, 75.1/scale], [-0.1/scale, 150/scale], color=wall_color, linewidth=wall_linewidth, alpha=0.7)

    if add_director:
        q_name = f"{save_dir}/"
        q_name += f"{file.replace('/','_')}_Q.npy"
        addNematicDirector(ax, cells, q_name, streamplot=streamplot, dr=3)

    ax.axis('scaled')
    ax.axis('off')
    # # Adjust plot limits to fit the walls
    # ax.set_xlim([-8, 80])  # Add a bit of margin
    # ax.set_ylim([-0.5, 0.2])  # Add a bit of margin

    ax.set_xlim([minx, maxx])
    ax.set_ylim([miny-0.1, maxy+0.1])

    stem = f"{save_dir}/plotCells_{file.replace('/','_')}"
    fig.savefig(f"{stem}_{streamplot=}.png",
                dpi=300,
                transparent=False,
                bbox_inches='tight')
    fig.savefig(f"{stem}_{streamplot=}.pdf",
                transparent=False,
                bbox_inches='tight')
    # plt.show()

if __name__ == "__main__":
    if int(sys.argv[2]) == 0:
        add_director = False
    else:
        add_director = True
    plotCells(sys.argv[1], add_director=add_director)