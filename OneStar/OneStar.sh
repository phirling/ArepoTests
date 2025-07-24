#!/bin/bash            # this line only there to enable syntax highlighting in this file

##################################################
#  Enable/Disable compile-time options as needed #
##################################################

OVERRIDE_PEANOGRID_WARNING
GADGET2_HEADER                                   # allows Arepo to understand ancient header formats

#--------------------------------------- Basic operation mode of code
NTYPES=6                                 # number of particle types

#GENERIC_ASYNC                           # enables asynchronous communication scheme
#PERIODIC

#USE_DIRECT_IO_FOR_RESTARTS 

MHD
MHD_POWELL
MHD_SEEDFIELD
MHD_POWELL_LIMIT_TIMESTEP

#--------------------------------------- Mesh Type
VORONOI

#--------------------------------------- Riemann solver
RIEMANN_HLLD

#--------------------------------------- Mesh motion and regularization
REGULARIZE_MESH_CM_DRIFT
REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED
REGULARIZE_MESH_FACE_ANGLE

#--------------------------------------- Time integration options
TREE_BASED_TIMESTEPS         # non-local timestep criterion (take 'signal speed' into account)

#--------------------------------------- Refinement and derefinement
REFINEMENT_SPLIT_CELLS
REFINEMENT_MERGE_CELLS

#REFINEMENT_HIGH_RES_GAS

#--------------------------------------- Gravity softening
NSOFTTYPES=6                      # Number of different softening values to which particle types can be mapped.
MULTIPLE_NODE_SOFTENING           # If a tree node is to be used which is softened, this is done with the softenings of its different mass components
ADAPTIVE_HYDRO_SOFTENING

#--------------------------------------- Things that are always recommended
CHUNKING                     # will calculated the gravity force in interleaved blocks. This can reduce imbalances in case multiple iterations due to insufficient buffer size need to be done

#---------------------------------------- Single/Double Precision
DOUBLEPRECISION=1
DOUBLEPRECISION_FFTW
# IDS_64BIT # On in LtU branch < ------- !!!

OUTPUT_COORDINATES_IN_DOUBLEPRECISION # to be implemented
NGB_TREE_DOUBLEPRECISION  # if this is enabled, double precision is used for the neighbor node extension

#---------------------------------------- On the fly FOF groupfinder
#FOF                                    # enable FoF output
#FOF_PRIMARY_LINK_TYPES=2               # 2^type for the primary dark matter type
#FOF_SECONDARY_LINK_TYPES=1+16+32   # 2^type for the types linked to nearest primaries
#FOF_SECONDARY_LINK_TARGET_TYPES=2+4+8   # should normally be set to a list of all dark matter types (in zoom runs), if not set defaults to FOF_PRIMARY_LINK_TYPES
#FOF_GROUP_MIN_LEN=32                   # default is 32
#FOF_LINKLENGTH=0.16                    # Linkinglength for FoF (default=0.2)
#FOF_FUZZ_SORT_BY_NEAREST_GROUP=0   # sort fuzz particles by nearest group and generate offset table in catalog (=1 writes nearest group number to snapshot)
#FOF_STOREIDS                           # store IDs in group/subfind catalogue, do not order particles in snapshot files by group order
#USE_AREPO_FOF_WITH_GADGET_FIX          # Needed in order to run FOF with Arepo on Gadget snapshot files, if gas is present and should be linked to the FOFs
#ADD_GROUP_PROPERTIES                   # This can be used to calculate additional properties for an already existing group catalogue. These are then added as additional columns to the HDF5 group catalogues.
#ADD_SO_GROUP_PROPERTIES                # This can be used to calculate additional properties for an already existing group catalogue. These are then added as additional columns to the HDF5 group catalogues.


#---------------------------------------- Subfind
#SUBFIND                                # enables substructure finder
#SUBFIND_HBT
#MERGERTREE
#SUBFIND_MEASURE_H2MASS                 # special measurement option for mass in molecular hydrogen
#SUBFIND_CALC_MORE                      # calculates also the velocity dispersion in the local density estimate
#SUBFIND_EXTENDED_PROPERTIES            # adds calculation of further quantities related to angular momentum in different components

#-------------------------------------------- Things for special behaviour
#READ_DM_AS_GAS
#NO_ID_UNIQUE_CHECK
#RUNNING_SAFETY_FILE                    # if file './running' exists, do not start the run
#LOAD_TYPES=1+2+4+16+32
#READ_COORDINATES_IN_DOUBLE
#IDS_OFFSET=1           #offset for gas particles if created from DM
#TILE_ICS
#COMBINETYPES            # reads in the IC file types 4+5 as type 3 (useful for doing gas runs of Aquarius ICs)
#USE_RANDOM_GENERATOR
#MULTIPLE_RESTARTS
#TOLERATE_WRITE_ERROR
#OPTIMIZE_MEMORY_USAGE                       #optimize for memory, not for speed. Note: this is dangerous for high dynamic range simulations with mixed precision, since some position variables are singles instead of doubles
#SUBBOX_SNAPSHOTS # note: maybe one only subbox0 of Illustris, better time spacing, smaller number of chunks (to be made a new param separated from the standard snaps)
### PROCESS_TIMES_OF_OUTPUTLIST
#EXTENDED_GHOST_SEARCH               # This extends the ghost search to the full 3x3 domain instead of the principal domain
#DOUBLE_STENCIL                     # this will ensure that the boundary region of the local mesh is deep enough to have a valid double stencil for all local cells
#TETRA_INDEX_IN_FACE                # adds an index to each entry of VF[] and DC[] to one of the tetrahedra that share this edge
VORONOI_DYNAMIC_UPDATE              # keeps track of mesh connectivity, which speeds up mesh construction
#COFFEE_PROBLEM
#NOH_PROBLEM
#SHIFT_BY_HALF_BOX
#DISABLE_VELOCITY_CSND_SLOPE_LIMITING
NO_MPI_IN_PLACE # On in LtU branch < ------- !!!
NO_ISEND_IRECV_IN_DOMAIN # On in LtU branch < ------- !!!
FIX_PATHSCALE_MPI_STATUS_IGNORE_BUG # On in LtU branch < ------- !!!
#USE_MPIALLTOALLV_IN_DOMAINDECOMP
#MPI_HYPERCUBE_ALLGATHERV         # some MPI-libraries may use quite a bit of internal storage for MPI_Allgatherv. This uses hypercubes instead as a work-around
#MPISENDRECV_CHECKSUM
#NOTREERND
ENLARGE_DYNAMIC_RANGE_IN_TIME   # This extends the dynamic range of the integer timeline from 32 to 64 bit Likely for the smaller box IllustrisDwarf
### NOSTOP_WHEN_BELOW_MINTIMESTEP
#DO_NOT_CREATE_STAR_PARTICLES
#DMPIC                              # enable special image code for dark matter simulations   
#ALLOWEXTRAPARAMS
#RADIATIVE_RATES                   # used in non-equilibrium chemistry model
#FIX_SPH_PARTICLES_AT_IDENTICAL_COORDINATES  # this can be used to load SPH ICs that contain identical particle coordinates
#VEL_POWERSPEC                     # compiles in a code module that allows via restart-flag 7 the calculation of a gas velocity power spectrum of a snapshot
#ADJ_BOX_POWERSPEC             # compiles in a code module that allows via restart-flag 7 the calculation of gas power spectra of a snapshot with an adjustable box (user defined center and size)
#DISABLE_OPTIMIZE_DOMAIN_MAPPING

#--------------------------------------- Output/Input options
#UPDATE_GRADIENTS_FOR_OUTPUT
REDUCE_FLUSH #commented for Eiger
#OUTPUT_REFBHCOUNTER                 
#OUTPUT_EVERY_STEP
#GODUNOV_STATS
OUTPUT_CPU_CSV
#OUTPUT_TASK
#OUTPUT_TIMEBIN_HYDRO
#OUTPUT_PRESSURE_GRADIENT
#OUTPUT_DENSITY_GRADIENT
#OUTPUT_VELOCITY_GRADIENT
#OUTPUT_BFIELD_GRADIENT
#OUTPUT_VERTEX_VELOCITY
#OUTPUT_VERTEX_VELOCITY_DIVERGENCE  # requires CALCULATE_VERTEX_VELOCITY_DIVERGENCE
OUTPUT_CENTER_OF_MASS
#OUTPUT_SURFACE_AREA
#OUTPUT_PRESSURE
#OUTPUTPOTENTIAL
#OUTPUTACCELERATION
#OUTPUTTIMESTEP
#OUTPUT_SOFTENINGS                # output particle softenings
#OUTPUTGRAVINTERACTIONS           # output gravitatational interactions (from the tree) of particles
HAVE_HDF5                         # needed when HDF5 I/O support is desired
#PARAMS_IN_SNAP                    # add the compiler flags and parameter file values to every snapshot file (requires HAVE_HDF5)
#HDF5_FILTERS                      # activate snapshot compression and checksum for HDF5 output
#OUTPUT_XDMF                       #writes an .xmf file for each snapshot, which can be read by visit (with the hdf5 snapshot)
#OUTPUTCOOLRATE                    # outputs cooling rate, and conduction rate if enabled
#OUTPUT_DIVVEL                             # output  velocity divergence
#OUTPUT_CURLVEL                     # output  velocity curl
#OUTPUT_COOLHEAT                   # output actual energy loss/gain in cooling/heating routine
#OUTPUT_VORTICITY
#OUTPUT_CELL_SPIN                  
#MEASURE_DISSIPATION_RATE          # measures and outputs dissipation rate. Note: requires USE_ENTROPY_FOR_COLD_FLOWS, even though it will then always use the thermal energy update
#OUTPUT_MACHNUM                    # output maximum mach number of a cell

#--------------------------------------- Testing and Debugging options
DEBUG                             # enables core-dumps
#DEBUG_ENABLE_FPU_EXCEPTIONS   # tries to enable FPU exceptions
#CHECKSUM_DEBUG
#RESTART_DEBUG
#VERBOSE                           # reports readjustments of buffer sizes
HOST_MEMORY_REPORTING             # reports after start-up the available system memory by analyzing /proc/meminfo
#FORCETEST=0.001                   # calculates for given fraction of particles direct summation forces to check accuracy of tree force
#FORCETEST_TESTFORCELAW=1          # this enables a special test to measure the effective force law of the code, can be set to 1 or 2

#--------------------------------------- Glass making/ 2nd-order initial conditions / Initial conditions options
#SECOND_ORDER_ICS
OFFSET_FOR_NON_CONTIGUOUS_IDS      # This needs to be set, because in Arepo2 the ID default offset is otherwise too low, and we get same ID errors
#GENERATE_GAS_IN_ICS
#SPLIT_PARTICLE_TYPE=2+4+8  #2
#NTYPES_ICS=6 # number of particle types in ICs, if not NTYPES (only works for 6, and non-HDF5 ICs!)
#SKIPOMEGACHECK       #avoids checking cosmological parameters? added by Matteo

#--------------------------------------- SGChem chemistry module
SGCHEM
CHEMISTRYNETWORK=10
#MCMA
#ABHE
#CHEM_IMAGE
#IMAGE_FOOTERS
SGCHEM_VARIABLE_Z                     #Allow metallicity and dust-to-gas ratio to vary between different cells
#SGCHEM_VARIABLE_ISRF                  #Allow interstellar radiation field strength to vary spatially
#SGCHEM_VARIABLE_CRION                 #Allow cosmic ray ionization rate to vary spatially
SGCHEM_TEMPERATURE_FLOOR
#SGCHEM_ACCRETION_LUMINOSITY
#SGCHEM_NO_HIGHN_DCHEM
#SGCHEM_DUMP_THERMAL_RATES
#SGCHEM_NO_MOLECULES
#SGCHEM_NO_COOL
#SGCHEM_SWITCHOFF_CHEM_ABOVE_THRESHOLD
SGCHEM_ADD_CONSTANT_PHOTORATES       #Specify spatially constant photoionisation and photoheating rates in addition to SWEEP-tracked values
CHEMCOOL
#NO_CO_COOL
#ONLY_DUST_COOL
#NO_DUST_H2_HEAT
#NO_GAS_DUST_HEAT
#ADIABATIC_DENSITY_THRESHOLD
#_13CO_
#COLLECT_EQB_DATA
#OLD_CO_SELF_SHIELDING
#THERMAL_INFO_DUMP
#DEBUG_SGCHEM
#DEBUG_COOLING_ID
#DEBUG_PARTICLE_ID
#DEBUGDUSTID
#DEBUG_RATE_EQ
#DEBUG_EVOLVE

#--------------------------------------- ECOGAL star formation prescription
SF_ECOGAL
SF_ECOGAL_CHECKS
SF_ECOGAL_LOG
#SF_ECOGAL_FORCE_SF_IN_DENS_GAS
# SF_ECOGAL_ADDITIONAL_OUTPUT
#SF_ECOGAL_STELLAR_TARGET_MASS
SF_ECOGAL_FEEDBACK
SF_ECOGAL_FEEDBACK_REINJECT_MASS
#SF_ECOGAL_FEEDBACK_PHOTOION
#SF_ECOGAL_POTENTIALPEAK
#SF_ECOGAL_FEEDBACK_PHOTOION_STROMGRENAPPROX
#SF_ECOGAL_FEEDBACK_PHOTOION_DEBUG
#REFINEMENT_AROUND_STARP_ECOGAL
#WARMUP_SN
MAXNUMSNE=1000
#PAT_DEBUG
#PAT_FORCE_MOMENTUM_INJ
PAT_SKIP_OMEGA_CHECK


#--------------------------------------- Supernova Feedback
SNE_FEEDBACK
#CLUSTERED_SNE
#INJECT_TRACER_INTO_SN
#SNE_solarring

#SNE_FEEDBACK_RANDOM_INCLUDE_STARP
#SNE_FEEDBACK_REDUCE_STELLAR_LIFETIMES
#SNE_FEEDBACK_REDUCE_FRACTION_HIGHMASS_STARS

PAT_RESAMPLE_IMF
PAT_DEBUG
#SF_ECOGAL_DISCRETE_ENRICHMENT
VERBOSE