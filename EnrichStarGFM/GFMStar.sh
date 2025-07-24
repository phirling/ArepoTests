#!/bin/bash            # this line only there to enable syntax highlighting in this file

##################################################
#  Enable/Disable compile-time options as needed #
##################################################

OVERRIDE_PEANOGRID_WARNING
GADGET2_HEADER                                   # allows Arepo to understand ancient header formats
PAT_DEBUG

#--------------------------------------- Basic operation mode of code
NTYPES=6                                 # number of particle types

#GENERIC_ASYNC                           # enables asynchronous communication scheme
#PERIODIC

#USE_DIRECT_IO_FOR_RESTARTS 

MHD
MHD_POWELL
MHD_SEEDFIELD
MHD_POWELL_LIMIT_TIMESTEP

COOLING
UVB_SELF_SHIELDING                # gas is self-shielded from the cosmic background based on its density
USE_SFR

PAT_SKIP_OMEGA_CHECK
VERBOSE

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

SOFTEREQS


GFM                                                        #master switch
GFM_STELLAR_EVOLUTION=0                    #stellar evolution: 0->default, 1->no mass loss (beta value changes + MassMetallicity & MassMetals inconsistent internally with cell dynamical mass) 2->call only test routine
GFM_CONST_IMF=0                            #0 for Chabrier (default), 1 for a pure power-law (requires parameter IMFslope, e.g. -2.35 for Salpeter)
#GFM_VARIABLE_IMF=0                         #0 for a pure power-law that depends on DM-veldisp
GFM_PREENRICH                              #pre enrich gas at given redshift
#GFM_EXACT_NUMNGB                           #use direct neighbor count instead of kernel weighted neighbor count
GFM_COOLING_METAL                          #metal line cooling
GFM_OUTPUT_MASK=1+2+4+8+16+32+64+256   #which fields to output (search GFM_OUTPUT_MASK in io.c to see which fields the bits encode)
#GFM_DUST                                   #formation and evolution of dust, requires GFM_STELLAR_EVOLUTION
#GFM_DUST_DESTMODE=0                        #dust destruction mode: 0->default (uses supernova rate), 1->constant destruction timescale
#GFM_CHECKS                                 #this checks the consistency of the AuxDataID/PID indices of stars and black holes every timestep
#GFM_DISCARD_ENRICHMENT_GRADIENTS           #this disables the gradient extrapolation of the passively advected metallicity scalar variables
GFM_NORMALIZED_METAL_ADVECTION             #this introduces an additional pseudo element for all untracked metals and normalizes the extrapolated abundance vectors to unity
GFM_OUTPUT_BIRTH_POS
GFM_AGN_RADIATION                          #cooling suppression/heating due to AGN radiation field (proximity effect)

#GFM_CHEMTAGS
#GFM_SPLITFE
#GFM_RPROCESS

GFM_DISCRETE_ENRICHMENT
#GFM_NO_NEGATIVE_ELEMENT_MASS_RELEASED  #do not allow that negative yields for each element consume more mass than contained in ejecta elemental composition


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
