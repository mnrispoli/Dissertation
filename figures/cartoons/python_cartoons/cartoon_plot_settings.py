__author__ = 'sine2k'

import pylab as pl
import numpy as np
#import brewer2mpl
from matplotlib.colors import LinearSegmentedColormap

#-------------------colormaps-------------------

#rawpic_cmap = brewer2mpl.get_map('RdBu', 'Diverging', 9).mpl_colormap

cmap = pl.get_cmap('RdBu')
colors = cmap(np.linspace(0.5, 1, cmap.N // 2))
density_colorbar = LinearSegmentedColormap.from_list('Upper Half', colors)

linestyles = { 
        "-":    (None, None), 
        "--":   (3., 2.), 
        "--.":   (3., 2., 3., 2., 1., 2.,),
        "---.":   (3., 2., 3., 2., 3., 2., 1., 2.,),
        "----.":   (3., 2., 3., 2., 3., 2., 3., 2., 1., 2.,),
        "-.":   (3., 2., 1., 2.,),
        "-..":   (3., 2., 1., 2., 1., 2.,),
        "-...":   (3., 2., 1., 2., 1., 2., 1., 2.),
        "-....":   (3., 2., 1., 2., 1., 2., 1., 2., 1., 2.),
        ":":    (1., 2.),
        }

#-------------------colors-------------------

# colors
color_winered = '#a40000'
color_red_m = "#BA0000"
color_red_d = '#660821'#"#720000"
# color_red_m = "#EE2B2D"
# color_red_d = "#A31D21"
color_red_l = '#B38190'#'#9C4C4C'

color2 = '#5c3566'
color_grey = "#7A7A7A"
color_black = '#000000'

color_violet_d = '#5c3566'
color_violet_m = '#75507b'
color_violet_l = '#cd7fa8'

color_blue_d = '#174074'#'#204a87'
color_blue_m = '#3465a4'
color_blue_l = '#859DBA'#'#729fef'

color_green_d = "#4e9a06"
color_green_m = "#73d216"

#color_brown_l = "#e9b96e"
color_brown_l = "#c4a000"
color_brown_m = "#774100"
color_brown_d = "#512E00"

dark_grey = "#4C4C4C"

#-------------------rcparams-------------------

rcpars = {
        'xtick.labelsize'      : 6,
        'xtick.major.size'     : 2,      # major tick size in points
        'xtick.minor.size'     : 1,      # minor tick size in points
        'xtick.major.width'    : 0.8,
        'xtick.minor.width'    : 0.4,
    
        'ytick.labelsize'      : 6,
        'ytick.major.size'     : 2,      # major tick size in points
        'ytick.minor.size'     : 1,      # minor tick size in points
        'ytick.major.width'    : 0.8,    # 0.5 pt minimum size for prl
        'ytick.minor.width'    : 0.4,
        
        'grid.linewidth'   :   0.5,

        'font.family' : 'serif',
        'font.sans-serif': 'Helvetica'
        
#        'text.usetex':True
#        
#        'font.family':'serif'
        
                
#        'xtick.labelsize': 6,
#        'ytick.labelsize': 6
                
        }

#rcpars = {
#
#    'font.sans-serif': 'Helvetica', #'Helvetica Neue'
#    'font.size': 6,
#    'mathtext.fontset': 'custom', #'dejavuserif', 'cm', 'custom', 'stix', 'stixsans', 'dejavusans'
#    'mathtext.rm': 'Helvetica Neue',
#    'mathtext.it': 'Helvetica Neue:italic',
#    'mathtext.bf': 'Helvetica Neue:bold',
#
#    # "font.serif": ["Helvetica Neue", "Bitstream Vera Sans"],
#    #"font.serif": ["Helvetica Neue", "Bitstream Vera Sans"],
#    #"font.serif": "Bembo Std",
#
#    # "font.sans-serif": ["Helvetica", "URW Nimbus Sans", "Bitstream Vera Sans"],
#    "font.family": "serif",  #"sans-serif",
#    "font.monospace": ["Courier", "Bitstream Vera Sans"],
#    #"font.family": "Gentium",
#
#    "pdf.fonttype": 42,  # correct cid subset, but problems with
#                              # glyphs
#    #"text.usetex": True, # uses better type1 fonts but but blows up
#                              # file size by 10
#
#    # activate this to use sans serif fonts in math mode in combination with text.usetex=true
##    "text.latex.preamble": [r"\usepackage{sfmath}"],
##    "text.usetex": True,
#
#    #"mathtext.default": "regular",
#
#    # "mathtext.fontset": "custom",
#    # "mathtext.cal": "cursive",
#    # "mathtext.rm": "serif",
#    # "mathtext.tt": "monospace",
#    # "mathtext.it": "serif:oblique", #"serif:italic",
#    # "mathtext.bf": "serif:bold",
#    # "mathtext.sf": "serif",
#    # "mathtext.fallback_to_cm": True,
#
#    "patch.linewidth": 0.5,
#
#    "figure.subplot.left": .25,
#    "figure.subplot.right": .94,
#    "figure.subplot.bottom": .2,
#    "figure.subplot.top": .95,
#    "figure.subplot.wspace": .2,
#    "figure.subplot.hspace": .2,
#
#
#    #'backend': 'ps',
#    'pdf.fonttype': 42,
#
#    'axes.labelsize': 6,
#    # 'axes.elinewidth': 0.5,
#    'axes.linewidth' : 0.4,
#
#    'lines.markersize': 2.5,
#    'lines.linewidth': 0.8,
#
#    'legend.fontsize': 5,
#    'legend.frameon': False,
#    'xtick.labelsize': 6,
#    'ytick.labelsize': 6,
#
#    "grid.linewidth"   :   0.5,
#
#    "xtick.major.size"     : 1.5,      # major tick size in points
#    "xtick.minor.size"     : 0.8,      # minor tick size in points
#    'xtick.major.width'    : 0.3,
#    'xtick.minor.width'    : 0.2,
#
#    "ytick.major.size"     : 1.5,      # major tick size in points
#    "ytick.minor.size"     : 0.8,      # minor tick size in points
#    'ytick.major.width'    : 0.3,
#    'ytick.minor.width'    : 0.2,
#
##
##
## 	"font.size": 8,
##
## 	"font.family": "cmr10",
## #	"font.family": "Gentium",
## 	"pdf.fonttype": 42,  # correct cid subset, but problems with
##                               # glyphs
## #	"text.usetex": True, # uses better type1 fonts but but blows up
##                               # file size by 10
## 	"mathtext.fontset": "cm",
## #	"mathtext.cal": "cursive",
## #	"mathtext.rm": "serif",
## #	"mathtext.tt": "monospace",
## #	"mathtext.it": "serif:italic",
## #	"mathtext.bf": "serif:bold",
## #	"mathtext.sf": "sans",
## 	"mathtext.fallback_to_cm": True,
##     "patch.linewidth": 0.5,
## 	"lines.linewidth": 0.5,
##     "lines.markeredgewidth": 0.5,
##     "lines.markersize": 3,
## 	"axes.labelsize": 8,
## 	"axes.titlesize": 8,
## 	"axes.linewidth": 0.5,
## 	"text.fontsize": 8,
## 	"legend.fontsize": 8,
## 	"xtick.labelsize": 8,
## 	"ytick.labelsize": 8,
## 	"figure.subplot.left": .08,
## 	"figure.subplot.right": .98,
## 	"figure.subplot.bottom": .1,
## 	"figure.subplot.top": .98,
## 	"figure.subplot.wspace": .2,
## 	"figure.subplot.hspace": .2,
## 	"legend.borderpad": .2,
## 	"legend.handlelength": 2.5,
## 	"legend.handletextpad": .01,
## 	"legend.borderaxespad": .3,
## 	# "legend.labelsep": .002,
## 	"legend.labelspacing": .1,
##     "axes3d.grid": False,
#}