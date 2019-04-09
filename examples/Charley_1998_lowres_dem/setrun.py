# encoding: utf-8
"""
Module to set up run time parameters for Clawpack.

The values set in the function setrun are then written out to data files
that will be read in by the Fortran code.

"""

from __future__ import absolute_import
from __future__ import print_function

import os
import datetime
import shutil
import gzip

import numpy as np
from geopy.distance import geodesic

from clawpack.geoclaw.surge.storm import Storm
import clawpack.clawutil as clawutil


def licricize_rmw(t, storm):
    """Eqn from p.3 of LICRICE docs"""
    # if RMW recorded at some timesteps
    # interpolate
    if storm.max_wind_radius.max() > 0:
        sr_valid = storm.max_wind_radius[storm.max_wind_radius>0]
        t_valid = np.array(storm.t)[storm.max_wind_radius>0]
        t_sec = [(s - datetime.datetime(2000,1,1)).total_seconds() for s in t_valid]
        this_t = (t - datetime.datetime(2000,1,1)).total_seconds()
        return np.interp(this_t,t_sec,sr_valid)

    # otherwise, if all RMW missing, guess, except
    # if there's no central pressure measurement
    ix = storm.t.index(t)
    if (np.nanmax(storm.central_pressure) > 0) and ((storm.central_pressure[ix] < 0) or (np.isnan(storm.central_pressure[ix]))):
        return -1

    vmax_obs = storm.max_wind_speed[ix]
    lat = storm.eye_location[ix,1]

    # remove translational speed
    # use centered difference to calc
    # translational speed if not at beginning
    # or end of storm
    if ix == 0:
        tsteps = [ix]
    elif ix == len(storm.t)-1:
        tsteps = [ix-1]
    else:
        tsteps = [ix-1,ix]
    v_trans = []
    for i in tsteps:
        this_x = storm.eye_location[i]
        next_x = storm.eye_location[i+1]
        # flip lon/lat to lat/lon
        dx = geodesic((this_x[1],this_x[0]),(next_x[1],next_x[0])).meters

        dt = (storm.t[i+1] - storm.t[i]).total_seconds()
        v_trans.append(dx/dt)
    v_trans = np.array(v_trans).mean()
    vmax_az = vmax_obs - v_trans

    return 63273. - 868.3*vmax_az + 1070. * lat


def set_missing_storm_rad_to_1000(t, storm):
        """
        Sets the default maximum storm radius to 1000m

        Holland 1980 parameterization does not decay to 0, so
        the geoclaw implementation (in
        ``geoclaw/src/2d/shallow/surge/model_storm_module.f90``)
        constrains windspeeds to being within 2 * the storm radius
        """
        return 1000


TRACK_FORMAT = 'IBTrACS'
TRACK_PATH = '/home/examples/data/IBTrACS.NA.hotel1.nc'
STORM_KWARGS = {'storm_name': 'CHARLEY', 'year': 1998}
RECURRENCE = 4
SL_INIT = -0.04563426870748305
FIXEDGRIDS = []
FGMAX_FILES = []
CHECKPOINT_STYLE = 0
TOPOFILES = [
    [3, 1, 5, -57600.0, 136800.0, '/home/examples/data/charlie_1998_base_dem.tty']]
LOWEST_RES = .25 # degree
T0 = -57600.0
TFINAL = 136800.0
DOMAIN = (-106.2, 17.75, -59.95, 46.25)
GAUGES = [[0, -94.985, 29.6817, T0, TFINAL], [1, -95.2658, 29.7262, T0, TFINAL], [2, -94.9183, 29.48, T0, TFINAL], [3, -94.7933, 29.31, T0, TFINAL], [4, -97.0467, 28.0217, T0, TFINAL], [5, -97.2164, 27.5808, T0, TFINAL], [6, -97.2155, 26.0612, T0, TFINAL]]
REGIONS = []
MIN_GAUGE_TIME_INCREMENT = 3600.0
REFINEMENT_RATIOS = [2, 2, 2, 6]
CHECKPOINT_INTERVAL = 1
CHECKPOINT_TIMES = []
RESTART = False
RS_FILE = ''
VERBOSITY = 0
RMW_FUNC = licricize_rmw
STORM_RAD_FUNC = set_missing_storm_rad_to_1000
PLOTS = ['surface', 'friction', 'pressure', 'gauge_vals', 'gauge_locs', 'wind', 'speed']


# Scratch directory for storing topo and dtopo files:
scratch_dir = os.path.join(os.environ["CLAW"], 'geoclaw', 'scratch')


# ------------------------------
def setrun(claw_pkg='geoclaw'):

    """
    Define the parameters used for running Clawpack.

    INPUT:
        claw_pkg expected to be "geoclaw" for this setrun.

    OUTPUT:
        rundata - object of class ClawRunData

    """

    from clawpack.clawutil import data

    assert claw_pkg.lower() == 'geoclaw',  "Expected claw_pkg = 'geoclaw'"

    num_dim = 2
    rundata = data.ClawRunData(claw_pkg, num_dim)

    # ------------------------------------------------------------------
    # Standard Clawpack parameters to be written to claw.data:
    #   (or to amr2ez.data for AMR)
    # ------------------------------------------------------------------
    clawdata = rundata.clawdata  # initialized when rundata instantiated

    # Set single grid parameters first.
    # See below for AMR parameters.

    # ---------------
    # Spatial domain:
    # ---------------

    # Number of space dimensions:
    clawdata.num_dim = num_dim

    # Lower and upper edge of computational domain:
    clawdata.lower[0] = DOMAIN[0]      # west longitude
    clawdata.upper[0] = DOMAIN[2]      # east longitude

    clawdata.lower[1] = DOMAIN[1]       # south latitude
    clawdata.upper[1] = DOMAIN[3]      # north latitude

    # Number of grid cells:
    clawdata.num_cells[0] = int((clawdata.upper[0] - clawdata.lower[0]) /
        LOWEST_RES)
    clawdata.num_cells[1] = int((clawdata.upper[1] - clawdata.lower[1]) /
        LOWEST_RES)

    # ---------------
    # Size of system:
    # ---------------

    # Number of equations in the system:
    clawdata.num_eqn = 3

    # Number of auxiliary variables in the aux array (initialized in setaux)
    # First three are from shallow GeoClaw, fourth is friction and last 3 are
    # storm fields
    clawdata.num_aux = 3 + 1 + 3

    # Index of aux array corresponding to capacity function, if there is one:
    clawdata.capa_index = 2

    # -------------
    # Initial time:
    # -------------
    clawdata.t0 = T0

    # Restart from checkpoint file of a previous run?
    # If restarting, t0 above should be from original run, and the
    # restart_file 'fort.chkNNNNN' specified below should be in
    # the OUTDIR indicated in Makefile.

    clawdata.restart = RESTART               # True to restart from prior results
    clawdata.restart_file = RS_FILE  # File to use for restart data

    # -------------
    # Output times:
    # --------------

    # Specify at what times the results should be written to fort.q files.
    # Note that the time integration stops after the final output time.
    # The solution at initial time t0 is always written in addition.

    clawdata.output_style = 1
    # Output nout frames at equally spaced times up to tfinal:
    clawdata.tfinal = TFINAL
    recurrence = RECURRENCE
    clawdata.num_output_times = int((clawdata.tfinal - clawdata.t0) *
                                        recurrence / (60**2 * 24))
    if clawdata.restart:
        clawdata.output_t0 = True # output at initial (or restart) time?
    else:
        clawdata.output_t0 = True

    clawdata.output_format = 'binary'      # 'ascii' or 'netcdf'
    clawdata.output_q_components = 'none'   # could be list such as [True,True]
    if len(PLOTS) > 1:
        # 'surface' gets output in q, but all other relevant
        # plottable variables are in aux, so if they are in PLOTS
        # make sure to output all aux components
        clawdata.output_aux_components = 'all'
    else:
        clawdata.output_aux_components = 'none'
    clawdata.output_aux_onlyonce = False    # output aux arrays only at t0

    # ---------------------------------------------------
    # Verbosity of messages to screen during integration:
    # ---------------------------------------------------

    # The current t, dt, and cfl will be printed every time step
    # at AMR levels <= verbosity.  Set verbosity = 0 for no printing.
    #   (E.g. verbosity == 2 means print only on levels 1 and 2.)
    clawdata.verbosity = VERBOSITY

    # --------------
    # Time stepping:
    # --------------

    # if dt_variable==1: variable time steps used based on cfl_desired,
    # if dt_variable==0: fixed time steps dt = dt_initial will always be used.
    clawdata.dt_variable = True

    # Initial time step for variable dt.
    # If dt_variable==0 then dt=dt_initial for all steps:
    clawdata.dt_initial = 0.016

    # Max time step to be allowed if variable dt used:
    clawdata.dt_max = 1e+99

    # Desired Courant number if variable dt used, and max to allow without
    # retaking step with a smaller dt:
    clawdata.cfl_desired = 0.75
    clawdata.cfl_max = 1.0

    # Maximum number of time steps to allow between output times:
    clawdata.steps_max = 2**16

    # ------------------
    # Method to be used:
    # ------------------

    # Order of accuracy:  1 => Godunov,  2 => Lax-Wendroff plus limiters
    clawdata.order = 1

    # Use dimensional splitting? (not yet available for AMR)
    clawdata.dimensional_split = 'unsplit'

    # For unsplit method, transverse_waves can be
    #  0 or 'none'      ==> donor cell (only normal solver used)
    #  1 or 'increment' ==> corner transport of waves
    #  2 or 'all'       ==> corner transport of 2nd order corrections too
    clawdata.transverse_waves = 1

    # Number of waves in the Riemann solution:
    clawdata.num_waves = 3

    # List of limiters to use for each wave family:
    # Required:  len(limiter) == num_waves
    # Some options:
    #   0 or 'none'     ==> no limiter (Lax-Wendroff)
    #   1 or 'minmod'   ==> minmod
    #   2 or 'superbee' ==> superbee
    #   3 or 'mc'       ==> MC limiter
    #   4 or 'vanleer'  ==> van Leer
    clawdata.limiter = ['mc', 'mc', 'mc']

    clawdata.use_fwaves = True    # True ==> use f-wave version of algorithms

    # Source terms splitting:
    #   src_split == 0 or 'none'
    #      ==> no source term (src routine never called)
    #   src_split == 1 or 'godunov'
    #      ==> Godunov (1st order) splitting used,
    #   src_split == 2 or 'strang'
    #      ==> Strang (2nd order) splitting used,  not recommended.
    clawdata.source_split = 'godunov'

    # --------------------
    # Boundary conditions:
    # --------------------

    # Number of ghost cells (usually 2)
    clawdata.num_ghost = 2

    # Choice of BCs at xlower and xupper:
    #   0 => user specified (must modify bcN.f to use this option)
    #   1 => extrapolation (non-reflecting outflow)
    #   2 => periodic (must specify this at both boundaries)
    #   3 => solid wall for systems where q(2) is normal velocity

    clawdata.bc_lower[0] = 'extrap'
    clawdata.bc_upper[0] = 'extrap'

    clawdata.bc_lower[1] = 'extrap'
    clawdata.bc_upper[1] = 'extrap'

    # Specify when checkpoint files should be created that can be
    # used to restart a computation.

    clawdata.checkpt_style = CHECKPOINT_STYLE

    if clawdata.checkpt_style == 2:
        clawdata.checkpt_times = CHECKPOINT_TIMES
    elif clawdata.checkpt_style == 3:
        clawdata.checkpt_interval = CHECKPOINT_INTERVAL

    # ---------------
    # AMR parameters:
    # ---------------
    amrdata = rundata.amrdata

    # max number of refinement levels:
    amrdata.amr_levels_max = len(REFINEMENT_RATIOS) + 1

    # List of refinement ratios at each level (length at least mxnest-1)
    amrdata.refinement_ratios_x = REFINEMENT_RATIOS
    amrdata.refinement_ratios_y = REFINEMENT_RATIOS
    amrdata.refinement_ratios_t = REFINEMENT_RATIOS

    # Specify type of each aux variable in amrdata.aux_type.
    # This must be a list of length maux, each element of which is one of:
    #   'center',  'capacity', 'xleft', or 'yleft'  (see documentation).

    amrdata.aux_type = ['center', 'capacity', 'yleft', 'center', 'center',
                        'center', 'center']

    # Flag using refinement routine flag2refine rather than richardson error
    amrdata.flag_richardson = False    # use Richardson?
    amrdata.flag2refine = True

    # steps to take on each level L between regriddings of level L+1:
    amrdata.regrid_interval = 3

    # width of buffer zone around flagged points:
    # (typically the same as regrid_interval so waves don't escape):
    amrdata.regrid_buffer_width = 3

    # clustering alg. cutoff for (# flagged pts) / (total # of cells refined)
    # (closer to 1.0 => more small grids may be needed to cover flagged cells)
    amrdata.clustering_cutoff = 0.700000

    # print info about each regridding up to this level:
    amrdata.verbosity_regrid = 0

    #  ----- For developers -----
    # Toggle debugging print statements:
    amrdata.dprint = False      # print domain flags
    amrdata.eprint = False      # print err est flags
    amrdata.edebug = False      # even more err est flags
    amrdata.gprint = False      # grid bisection/clustering
    amrdata.nprint = False      # proper nesting output
    amrdata.pprint = False      # proj. of tagged points
    amrdata.rprint = False      # print regridding summary
    amrdata.sprint = False      # space/memory output
    amrdata.tprint = False      # time step reporting each level
    amrdata.uprint = False      # update/upbnd reporting

    # More AMR parameters can be set -- see the defaults in pyclaw/data.py

    # == setregions.data values ==
    rundata.regiondata.regions = REGIONS
    rundata.gaugedata.gauges = GAUGES
    rundata.gaugedata.q_out_fields = [0]
    rundata.gaugedata.min_time_increment = MIN_GAUGE_TIME_INCREMENT


    # ------------------------------------------------------------------
    # GeoClaw specific parameters:
    # ------------------------------------------------------------------
    rundata = setgeo(rundata)

    return rundata
    # end of function setrun
    # ----------------------


# -------------------
def setgeo(rundata):
    """
    Set GeoClaw specific runtime parameters.
    For documentation see ....
    """

    geo_data = rundata.geo_data

    # == Physics ==
    geo_data.gravity = 9.81
    geo_data.coordinate_system = 2
    geo_data.earth_radius = 6367.5e3
    geo_data.rho = 1025.0
    geo_data.rho_air = 1.15
    geo_data.ambient_pressure = 101.3e3

    # == Forcing Options
    geo_data.coriolis_forcing = True
    geo_data.friction_forcing = True
    geo_data.friction_depth = 1e10

    # == Algorithm and Initial Conditions ==
    geo_data.sea_level = SL_INIT
    geo_data.dry_tolerance = 1.0e-2

    # Refinement Criteria
    refine_data = rundata.refinement_data
    refine_data.wave_tolerance = 1.0
    refine_data.speed_tolerance = [1.0, 2.0, 3.0, 4.0]
    refine_data.deep_depth = 300.0
    refine_data.max_level_deep = 4
    refine_data.variable_dt_refinement_ratios = True

    # == settopo.data values ==
    topo_data = rundata.topo_data
    topo_data.topofiles = TOPOFILES

    # == setfixedgrids.data values ==
    rundata.fixed_grid_data.fixedgrids = FIXEDGRIDS
    rundata.fgmax_data.fgmax_files = FGMAX_FILES
    rundata.fgmax_data.num_fgmax_val = 1

    # for fixed grids append lines of the form
    # [t1,t2,noutput,x1,x2,y1,y2,xpoints,ypoints,\
    #  ioutarrivaltimes,ioutsurfacemax]

    # ================
    #  Set Surge Data
    # ================
    data = rundata.surge_data

    # Source term controls
    data.wind_forcing = True
    data.drag_law = 1
    data.pressure_forcing = True

    data.display_landfall_time = True

    # AMR parameters, m/s and m respectively
    data.wind_refine = [20.0, 40.0, 60.0]
    data.R_refine = [60.0e3, 40e3, 20e3]

    # Storm parameters - Parameterized storm (Holland 1980)
    data.storm_specification_type = 1  # (type 1)
    data.storm_file = os.path.expandvars(os.path.join(os.getcwd(),
                                'storm.storm'))

    # Convert IBTrACS data to GeoClaw format
    storm = Storm(path=TRACK_PATH, file_format=TRACK_FORMAT,
               **STORM_KWARGS)
    storm.write(data.storm_file, file_format='geoclaw',
               max_wind_radius_fill = RMW_FUNC,
               storm_radius_fill = STORM_RAD_FUNC)

    # =======================
    #  Set Variable Friction
    # =======================
    data = rundata.friction_data

    # Variable friction
    data.variable_friction = True

    # Region based friction
    # Entire domain
    data.friction_regions.append([rundata.clawdata.lower,
                                  rundata.clawdata.upper,
                                  [np.infty, 0.0, -np.infty],
                                  [0.030, 0.025]])

    # La-Tex Shelf
    data.friction_regions.append([(-98, 25.25), (-90, 30),
                                  [np.infty, -10.0, -200.0, -np.infty],
                                  [0.030, 0.012, 0.022]])

    return rundata
    # end of function setgeo
    # ----------------------


if __name__ == '__main__':
    # Set up run-time parameters and write all data files.
    import sys
    if len(sys.argv) == 2:
        rundata = setrun(sys.argv[1])
    else:
        rundata = setrun()

    rundata.write()
