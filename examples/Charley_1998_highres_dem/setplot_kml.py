
"""
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.

"""


from __future__ import absolute_import
article = False

import os

import numpy

# Plot customization
import matplotlib

import matplotlib.pyplot as plt
import datetime

from clawpack.visclaw import colormaps, geoplot
import clawpack.clawutil.data as clawutil
import clawpack.amrclaw.data as amrclaw
import clawpack.geoclaw.data as geodata

import clawpack.geoclaw.surge.plot as surgeplot
import clawpack.geoclaw.surge.data as surgedata

try:
    from setplotfg import setplotfg
except:
    setplotfg = None
    
STORMNAME = 
TIME_OFFSET = '1998-08-22T10:00:00'
T0_DT = [1998, 8, 21, 18, 0, 0]
NOAA_GAUGE_NUMS = [0, 1, 2, 3, 4, 5, 6]
PROP_GAUGE_NUMS = []
CL_GAUGE_NUMS = []
GAUGE_NUMS = sum([NOAA_GAUGE_NUMS,PROP_GAUGE_NUMS,CL_GAUGE_NUMS],[])
INTERVAL_DAYS = [-1, 1]

def setplot(plotdata):
    r"""Setplot function for surge plotting"""


    plotdata.clearfigures()  # clear any old figures,axes,items data
    plotdata.format = 'binary'

    fig_num_counter = surgeplot.figure_counter()

    # Load data from output
    clawdata = clawutil.ClawInputData(2)
    clawdata.read(os.path.join(plotdata.outdir,'claw.data'))

    amrdata = amrclaw.AmrclawInputData(clawdata)
    amrdata.read(os.path.join(plotdata.outdir,'amr.data'))

    physics = geodata.GeoClawData()
    physics.read(os.path.join(plotdata.outdir,'geoclaw.data'))

    surge_data = surgedata.SurgeData()
    surge_data.read(os.path.join(plotdata.outdir,'surge.data'))

    friction_data = surgedata.FrictionData()
    friction_data.read(os.path.join(plotdata.outdir,'friction.data'))

    # Load storm track
    track = surgeplot.track_data(os.path.join(plotdata.outdir,'fort.track'))

    # Color limits
    surface_range = 5.0
    speed_range = 3.0
    eta = physics.sea_level
    if not isinstance(eta,list):
        eta = [eta]
    surface_limits = [eta[0]-surface_range,eta[0]+surface_range]
    # surface_contours = numpy.linspace(-surface_range, surface_range,11)
    surface_contours = [-5,-4.5,-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5]

    speed_limits = [0.0,speed_range]
    speed_contours = numpy.linspace(0.0,speed_range,13)

    wind_limits = [0,64]
    # wind_limits = [-0.002,0.002]
    pressure_limits = [935,1013]
    friction_bounds = [0.01,0.04]
    # vorticity_limits = [-1.e-2,1.e-2]

    # ==========================================================================
    # ==========================================================================
    #   Plot specifications
    # ==========================================================================
    # ==========================================================================

    #-----------------------------------------
    # Some global kml flags
    #-----------------------------------------
    plotdata.kml_name = STORMNAME
    plotdata.kml_starttime = T0_DT   # Time of event in UTC [None]
    #plotdata.kml_tz_offset = 0    # Time zone offset (in hours) of event. [None]

    plotdata.kml_index_fname = STORMNAME  # name for .kmz and .kml files ["_GoogleEarth"]

    # Set to path where KMZ files will be stored;   KML file will then
    # link to this path.
    # plotdata.kml_publish = 'http://www.domain.edu/path/to/kmz/files'

    # ========================================================================
    #  Entire Domain
    # ========================================================================
    full_xlimits = [clawdata.lower[0],clawdata.upper[0]]
    full_ylimits = [clawdata.lower[1],clawdata.upper[1]]
    full_shrink = 1.0

    # --------------------------
    #  Surface - entire domain
    # --------------------------

    plotfigure = plotdata.new_plotfigure(name='Surface - Entire Domain',
                                         figno=0)
    plotfigure.show = True

    plotfigure.use_for_kml = True
    plotfigure.kml_use_for_initial_view = True

    # These override axes limits set below in plotitems
    plotfigure.kml_xlimits = full_xlimits
    plotfigure.kml_ylimits = full_ylimits

    # Resolution - needs to be set carefully for the transparent colormap
    rcl = 40
    plotfigure.kml_dpi = rcl*2
    plotfigure.kml_figsize = [11.6, 9.6]
    plotfigure.kml_tile_images = False    # Tile images for faster loading.  Requires GDAL [False]

    plotaxes = plotfigure.new_plotaxes()
    plotitem = plotaxes.new_plotitem(name='surface',plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.surface_or_depth
    plotitem.pcolor_cmap = geoplot.googleearth_transparent
    plotitem.pcolor_cmin = -surface_range
    plotitem.pcolor_cmax = surface_range


    def kml_colorbar(filename):
        cmin = -surface_range
        cmax = surface_range
        geoplot.kml_build_colorbar(filename,
                                   geoplot.googleearth_transparent,
                                   cmin,cmax)

    plotfigure.kml_colorbar = kml_colorbar


    # --------------------------
    #  Water Speed - entire domain
    # --------------------------
    plotfigure = plotdata.new_plotfigure(name='Currents - Entire Domain',
                                         figno=1)
    plotfigure.show = True

    plotfigure.use_for_kml = True
    plotfigure.kml_use_for_initial_view = False

    # These override axes limits set below in plotitems
    plotfigure.kml_xlimits = full_xlimits
    plotfigure.kml_ylimits = full_ylimits
    plotfigure.kml_figsize = [11.6,9.6]
    plotfigure.kml_dpi = 80   # size not so important with contourf
    plotfigure.kml_tile_images = False    # Tile images for faster loading.  Requires GDAL [False]

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()

    plotitem = plotaxes.new_plotitem(name='surface', plot_type='2d_contourf')
    plotitem.plot_var = geoplot.surface_or_depth
    plotitem.contour_levels = speed_contours
    plotitem.fill_cmin = min(speed_contours)
    plotitem.fill_cmax = max(speed_contours)

    cmap= plt.get_cmap('OrRd')
    cmap._rgba_under = (0.0,0.0,0.0,0.0)
    plotitem.fill_cmap = cmap

    def cbar_speeds(filename):
        cmin = min(speed_contours)
        cmax = max(speed_contours)
        geoplot.kml_build_colorbar(filename,cmap,cmin,cmax)

    plotfigure.kml_colorbar = cbar_speeds


    # ========================================================================
    #  Figures for gauges
    # ========================================================================
    plotfigure = plotdata.new_plotfigure(name='Surface & topo', figno=300, \
                    type='each_gauge')
    plotfigure.show = True
    plotfigure.clf_each_gauge = True
    # plotfigure.kwargs['figsize'] = (16,10)

    def gauge_after_axes(cd):

        if cd.gaugeno in GAUGE_NUMS:
            axes = plt.gca()

            # Add GeoClaw gauge data
            geoclaw_gauge = cd.gaugesoln
            axes.plot([TIME_OFFSET + timedelta(seconds=t) for t in geoclaw_gauge.t],
                  geoclaw_gauge.q[1,:], 'b--',
                  label="GeoClaw")

            # Fix up plot
            axes.set_title('Station %s' % cd.gaugeno)
            axes.set_xlabel('Days relative to landfall')
            axes.set_ylabel('Surface (m)')
            axes.set_xlim(INTERVAL_DAYS)
            axes.set_ylim([-1,5])
            #axes.set_xticks([-2,-1,0,1])
            #axes.set_xticklabels([r"$-2$",r"$-1$",r"$0$",r"$1$"])
            axes.grid(True)
            axes.legend()

            plt.hold(False)

        surgeplot.gauge_afteraxes(cd)


    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = INTERVAL_DAYS
    plotaxes.ylimits = [-1,5]
    plotaxes.title = 'Surface'
    plotaxes.afteraxes = gauge_after_axes

    # Plot surface as blue curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 3
    plotitem.plotstyle = 'b-'

    #-----------------------------------------

    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_gaugenos = 'all'          # list of gauges to print
    plotdata.print_fignos = range(len(GAUGE_NUMS))            # list of figures to print
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?

    plotdata.html = False                     # create html files of plots?
    plotdata.latex = False                    # create latex file of plots?
    plotdata.kml = True

    return plotdata
