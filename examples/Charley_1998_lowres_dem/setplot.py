
from __future__ import absolute_import
from __future__ import print_function

import os

import numpy
import matplotlib.pyplot as plt
import datetime
import pandas as pd

import clawpack.visclaw.colormaps as colormap
import clawpack.visclaw.gaugetools as gaugetools
import clawpack.clawutil.data as clawutil
import clawpack.amrclaw.data as amrclaw
import clawpack.geoclaw.data as geodata


import clawpack.geoclaw.surge.plot as surgeplot

try:
    from setplotfg import setplotfg
except:
    setplotfg = None

NOAA_GAUGE_NUMS = [0, 1, 2, 3, 4, 5, 6]
PROP_GAUGE_NUMS = []
CL_GAUGE_NUMS = []
TIME_OFFSET = '1998-08-22T10:00:00'
PLOTS = ['surface', 'friction', 'pressure', 'gauge_vals', 'gauge_locs', 'wind', 'speed']
CREATE_PDF = False


def setplot(plotdata=None):
    """"""

    
    if plotdata is None:
        from clawpack.visclaw.data import ClawPlotData
        plotdata = ClawPlotData()

    # clear any old figures,axes,items data
    plotdata.clearfigures()
    plotdata.format = 'binary'

    # Load data from output
    clawdata = clawutil.ClawInputData(2)
    clawdata.read(os.path.join(plotdata.outdir, 'claw.data'))
    physics = geodata.GeoClawData()
    physics.read(os.path.join(plotdata.outdir, 'geoclaw.data'))
    surge_data = geodata.SurgeData()
    surge_data.read(os.path.join(plotdata.outdir, 'surge.data'))
    friction_data = geodata.FrictionData()
    friction_data.read(os.path.join(plotdata.outdir, 'friction.data'))

    # Load storm track
    track = surgeplot.track_data(os.path.join(plotdata.outdir, 'fort.track'))

    # Set afteraxes function
    def surge_afteraxes(cd):
        surgeplot.surge_afteraxes(cd, track, plot_direction=False,
                                             kwargs={"markersize": 4})

    # Color limits
    surface_limits = [-1.5,1.5]
    speed_limits = [0.0, 3.0]
    wind_limits = [0, 64]
    pressure_limits = [935, 1013]
    friction_bounds = [0.01, 0.04]

    def friction_after_axes(cd):
        plt.title(r"Manning's $n$ Coefficient")

    # ==========================================================================
    #   Plot specifications
    # ==========================================================================
    regions = {"Full Domain": {"xlimits": (clawdata.lower[0], clawdata.upper[0]),
                        "ylimits": (clawdata.lower[1], clawdata.upper[1]),
                        "figsize": (6.4, 4.8)}}

    for (name, region_dict) in regions.items():

        # Surface Figure
        if 'surface' in PLOTS:
            plotfigure = plotdata.new_plotfigure(name="Surface - %s" % name)
            plotfigure.kwargs = {"figsize": region_dict['figsize']}
            plotaxes = plotfigure.new_plotaxes()
            plotaxes.title = "Surface"
            plotaxes.xlimits = region_dict["xlimits"]
            plotaxes.ylimits = region_dict["ylimits"]
            plotaxes.afteraxes = surge_afteraxes

            surgeplot.add_surface_elevation(plotaxes, bounds=surface_limits)
            surgeplot.add_land(plotaxes)
            #plotaxes.plotitem_dict['surface'].amr_patchedges_show = [0] * 10
            #plotaxes.plotitem_dict['land'].amr_patchedges_show = [0] * 10

        # Speed Figure
        if 'speed' in PLOTS:
            plotfigure = plotdata.new_plotfigure(name="Currents - %s" % name)
            plotfigure.kwargs = {"figsize": region_dict['figsize']}
            plotaxes = plotfigure.new_plotaxes()
            plotaxes.title = "Currents"
            plotaxes.xlimits = region_dict["xlimits"]
            plotaxes.ylimits = region_dict["ylimits"]
            plotaxes.afteraxes = surge_afteraxes

            surgeplot.add_speed(plotaxes, bounds=speed_limits)
            surgeplot.add_land(plotaxes)
            #plotaxes.plotitem_dict['surface'].amr_patchedges_show = [0] * 10
            #plotaxes.plotitem_dict['land'].amr_patchedges_show = [0] * 10
        
    #
    # Friction field
    #
    if 'friction' in PLOTS:
        plotfigure = plotdata.new_plotfigure(name='Friction')
        plotfigure.show = friction_data.variable_friction and True

        plotaxes = plotfigure.new_plotaxes()
        plotaxes.xlimits = regions['Full Domain']['xlimits']
        plotaxes.ylimits = regions['Full Domain']['ylimits']
        # plotaxes.title = "Manning's N Coefficient"
        plotaxes.afteraxes = friction_after_axes
        plotaxes.scaled = True

        surgeplot.add_friction(plotaxes, bounds=friction_bounds, shrink=0.9)
        plotaxes.plotitem_dict['friction'].amr_patchedges_show = [0] * 10
        plotaxes.plotitem_dict['friction'].colorbar_label = "$n$"

    #
    #  Hurricane Forcing fields
    #
    # Pressure field
    if 'pressure' in PLOTS:
        plotfigure = plotdata.new_plotfigure(name='Pressure')
        plotfigure.show = surge_data.pressure_forcing and True

        plotaxes = plotfigure.new_plotaxes()
        plotaxes.xlimits = regions['Full Domain']['xlimits']
        plotaxes.ylimits = regions['Full Domain']['ylimits']
        plotaxes.title = "Pressure Field"
        plotaxes.afteraxes = surge_afteraxes
        plotaxes.scaled = True
        surgeplot.add_pressure(plotaxes, bounds=pressure_limits)
        surgeplot.add_land(plotaxes)

    # Wind field
    if 'wind' in PLOTS:
        plotfigure = plotdata.new_plotfigure(name='Wind Speed')
        plotfigure.show = surge_data.wind_forcing and True

        plotaxes = plotfigure.new_plotaxes()
        plotaxes.xlimits = regions['Full Domain']['xlimits']
        plotaxes.ylimits = regions['Full Domain']['ylimits']
        plotaxes.title = "Wind Field"
        plotaxes.afteraxes = surge_afteraxes
        plotaxes.scaled = True
        surgeplot.add_wind(plotaxes, bounds=wind_limits)
        surgeplot.add_land(plotaxes)

    # ========================================================================
    #  Figures for gauges
    # ========================================================================
    if 'gauge_vals' in PLOTS:
        plotfigure = plotdata.new_plotfigure(name='Gauge Surfaces', figno=300,
                                             type='each_gauge')
        plotfigure.show = True
        plotfigure.clf_each_gauge = True

        # Set up for axes in this figure:
        plotaxes = plotfigure.new_plotaxes()
        plotaxes.xlimits = 'auto'
        plotaxes.ylimits = 'auto'
        plotaxes.title = 'Surface'

        def gauge_afteraxes(cd):

            axes = plt.gca()
            fig = plt.gcf()
            gauge = cd.gaugesoln
            to = datetime.datetime.strptime(TIME_OFFSET,'%Y-%m-%dT%H:%M:%S')
            t = [to + datetime.timedelta(seconds=i) for i in gauge.t]
            axes.lines[0].set_xdata(t)
            axes.set_xlim([min(t),max(t)])
            fig.autofmt_xdate()
            fig.tight_layout()

            # Fix up plot - in particular fix time labels
            axes.set_title('Station %s' % cd.gaugeno)
            axes.set_xlabel('Datetime (UTC)')
            axes.set_ylabel('Surface (m)')
            axes.grid(True)
        plotaxes.afteraxes = gauge_afteraxes

        # Plot surface as blue curve:
        plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
        plotitem.plot_var = 1

    #
    #  Gauge Location Plot
    #
    if 'gauge_locs' in PLOTS:
        def gauge_location_afteraxes(cd):
            plt.subplots_adjust(left=0.12, bottom=0.06, right=0.97, top=0.97)
            surge_afteraxes(cd)
            
            # plot prop gauges
            gaugetools.plot_gauge_locations(cd.plotdata, gaugenos=PROP_GAUGE_NUMS,
                                            format_string='#1b9e77', add_labels=False,
                                            markersize=1)
            
            # plot coastline gauges
            gaugetools.plot_gauge_locations(cd.plotdata, gaugenos=CL_GAUGE_NUMS,
                                            format_string='#7570b3', add_labels=False,
                                            markersize=1)
            
            # plot noaa gauges
            gaugetools.plot_gauge_locations(cd.plotdata, gaugenos=NOAA_GAUGE_NUMS,
                                            format_string='#d95f02', add_labels=False,
                                            markersize=1)

        plotfigure = plotdata.new_plotfigure(name="Gauge Locations")
        plotfigure.show = True

        # Set up for axes in this figure:
        plotaxes = plotfigure.new_plotaxes()
        plotaxes.title = 'Gauge Locations'
        plotaxes.scaled = True
        plotaxes.xlimits = 'auto'
        plotaxes.ylimits = 'auto'
        plotaxes.afteraxes = gauge_location_afteraxes
        surgeplot.add_surface_elevation(plotaxes, bounds=surface_limits)
        surgeplot.add_land(plotaxes)
        plotaxes.plotitem_dict['surface'].amr_patchedges_show = [0] * 10
        plotaxes.plotitem_dict['land'].amr_patchedges_show = [0] * 10

    # -----------------------------------------
    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_gaugenos = NOAA_GAUGE_NUMS     # list of gauges to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.latex = CREATE_PDF              # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = CREATE_PDF      # also run pdflatex?
    plotdata.parallel = True                 # parallel plotting

    return plotdata