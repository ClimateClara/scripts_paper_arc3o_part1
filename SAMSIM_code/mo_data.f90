
!! Sets data and contains all flag descriptions.
!!
!! All data needed by mo_grotz are set in this module.
!! Most arrays are allocated after the needed dimension is specified for each testcase in mo_init.f90.
!!
!!
!! @author Philipp Griewank
!!
!!  COPYRIGHT
!!        Copyright (c) 2014 Max-Planck-Institut fuer Meteorologie, Hamburg,
!!        Germany
!!
!!        Copying and distribution of this file, with or without modification,
!!        are permitted in any medium without royalty provided the copyright
!!        notice and this notice are preserved.  This file is offered as-is,
!!        without any warranty.
!!
!!
!!
!!
!!
!! @par Revision History
!! Initialized by Philipp Griewank, IMPRS (2010-07-14)
!!

!!
MODULE mo_data

  USE mo_parameters, ONLY:wp

  IMPLICIT NONE

  
  !----Arrays
  REAL(wp), DIMENSION(:), ALLOCATABLE   :: H                 !<  Enthalpy [J]
  REAL(wp), DIMENSION(:), ALLOCATABLE   :: H_abs             !<  specific Enthalpy [J/kg]
  REAL(wp), DIMENSION(:), ALLOCATABLE   :: Q                 !<  Heat in layer [J]
  REAL(wp), DIMENSION(:), ALLOCATABLE   :: fl_Q              !<  Heat flux between layers [J/s]
  REAL(wp), DIMENSION(:), ALLOCATABLE   :: T                 !<  Temperature [C]
  REAL(wp), DIMENSION(:), ALLOCATABLE   :: S_bu              !<  Bulk Salinity [g/kg]
  REAL(wp), DIMENSION(:), ALLOCATABLE   :: fl_S              !<  Salinity flux [(g/s]
  REAL(wp), DIMENSION(:), ALLOCATABLE   :: S_abs             !<  Absolute Salinity [g]
  REAL(wp), DIMENSION(:), ALLOCATABLE   :: S_br              !<  Brine salinity [g/kg]
  REAL(wp), DIMENSION(:), ALLOCATABLE   :: thick             !<  Layer thickness [m]
  REAL(wp), DIMENSION(:), ALLOCATABLE   :: m                 !<  Mass [kg]
  REAL(wp), DIMENSION(:), ALLOCATABLE   :: fl_m              !<  Mass fluxes between layers [kg]
  REAL(wp), DIMENSION(:), ALLOCATABLE   :: V_s               !<  Volume [m^3] of solid
  REAL(wp), DIMENSION(:), ALLOCATABLE   :: V_l               !<  Volume [m^3] of liquid
  REAL(wp), DIMENSION(:), ALLOCATABLE   :: V_g               !<  Volume [m^3] of gas
  REAL(wp), DIMENSION(:), ALLOCATABLE   :: V_ex              !<  Volume of brine due expelled due to freezing [m^3] of solid, gas & liquid
  REAL(wp), DIMENSION(:), ALLOCATABLE   :: phi               !<  Solid mass fraction
  REAL(wp), DIMENSION(:), ALLOCATABLE   :: psi_s             !<  Solid volume fraction
  REAL(wp), DIMENSION(:), ALLOCATABLE   :: psi_l             !<  Liquid volume fraction
  REAL(wp), DIMENSION(:), ALLOCATABLE   :: psi_g             !<  Gas volume fraction
  REAL(wp), DIMENSION(:), ALLOCATABLE   :: ray               !<  Rayleigh number of each layer
  REAL(wp), DIMENSION(:), ALLOCATABLE   :: perm              !<  Permeability [?]

  REAL(wp)                              :: dt                !<  Timestep [s]
  REAL(wp)                              :: thick_0           !<  Initial layer thickness [m]
  REAL(wp)                              :: time              !<  Time [s]
  REAL(wp)                              :: freeboard         !<  Height of ice surface above (or below) waterlevel [m]
  REAL(wp)                              :: T_freeze          !<  Freezing temperature [C]
  INTEGER                               :: Nlayer            !<  Number of layers
  INTEGER                               :: N_bottom          !<  Number of bottom layers
  INTEGER                               :: N_middle          !<  Number of middle layers
  INTEGER                               :: N_top             !<  Number of top layers
  INTEGER                               :: N_active          !<  Number of Layers active in the present
  INTEGER                               :: i                 !<  Index, normally used for time
  INTEGER                               :: k                 !<  Index, normally used for layer
  REAL(wp)                              :: time_out          !<  Time between outputs [s]
  REAL(wp)                              :: time_total        !<  Time of simulation [s]
  INTEGER                               :: i_time            !<  Number of timesteps
  INTEGER                               :: i_time_out        !<  Number of timesteps between each output
  INTEGER                               :: n_time_out        !<  Counts number of timesteps between output
  CHARACTER*12000                       :: format_T,format_psi,format_thick,format_snow,format_integer,format_T2m_top,format_bgc  !< Format strings for output


  !----Boundary conditions
  REAL(wp)                              :: T_bottom          !<  Temperature of water beneath the ice [C]
  REAL(wp)                              :: T_top             !<  Temperature at the surface [C]
  REAL(wp)                              :: S_bu_bottom       !<  Salinity beneath the ice [g/kg]
  REAL(wp)                              :: T2m               !<  Two meter Temperature [C]
  REAL(wp)                              :: fl_q_bottom       !<  Bottom heat flux [J*s]

  !----Snow & Precip
  REAL(wp)                              :: psi_s_snow        !<  Solid volume fraction of snow layer
  REAL(wp)                              :: psi_l_snow        !<  Liquid volume fraction of snow layer
  REAL(wp)                              :: psi_g_snow        !<  Gas volume fraction of snow layer 
  REAL(wp)                              :: phi_s             !<  Solid mass fraction of snow layer 
  REAL(wp)                              :: S_abs_snow        !<  Absolute salinity of snow layer [g]
  REAL(wp)                              :: H_abs_snow        !<  Absolute enthalpy of snow layer [J]
  REAL(wp)                              :: m_snow            !<  Mass of snow layer [kg]
  REAL(wp)                              :: T_snow            !<  Temperature of snow layer [C]
  REAL(wp)                              :: thick_snow        !<  Thickness of snow layer [m]
  REAL(wp)                              :: liquid_precip     !<  Liquid precip, [meter of water/s]
  REAL(wp)                              :: solid_precip      !<  Solid precip, [meter of water /s]
  REAL(wp)                              :: fl_q_snow         !<  flow of heat into the snow layer  

  !----Vital signs
  REAL(wp)                              :: energy_stored     !<  Total amount of energy stored, control is freezing point temperature of S_bu_bottom [J]
  REAL(wp)                              :: total_resist      !<  Thermal resistance of the whole column []
  REAL(wp)                              :: surface_water     !<  Percentage of water fraction in the top 5cm [%]
  REAL(wp)                              :: freshwater        !<  Meters of freshwater stored in column [m]
  REAL(wp)                              :: thickness         !<  Meters of ice [m]
  REAL(wp)                              :: bulk_salin        !<  Salt/Mass [ppt]

  !----Model and numerics
  REAL(wp)                              :: thick_min         !<  Parameter for snow, determines when snow is in thermal equilibrium with the ice and when it is totally neglected
  REAL(wp),SAVE                         :: T_test            !<  First guess for getT subroutine

  !----Radiation
  REAL(wp)                              :: albedo            !<  Amount of short wave radiation which is reflected at the top surface
  REAL(wp)                              :: fl_sw             !<  Incoming shortwave radiation [W/m**2]
  REAL(wp)                              :: fl_lw             !<  Incoming longwave radiation  [W/m**2]
  REAL(wp)                              :: fl_sen            !<  Sensitive heat flux [W/m**2]
  REAL(wp)                              :: fl_lat            !<  Latent heat flux [W/m**2]
  REAL(wp)                              :: fl_rest           !<  Bundled longwave,sensitive and latent heat flux [W/m**2]
  REAL(wp), DIMENSION(:), ALLOCATABLE   :: fl_rad            !<  Energy flux of absorbed sw radiation of each layer [J/s]

  !----Gravity drainage
  REAL(wp)                              :: grav_drain        !<  brine flux of gravity drainage between two outputs [kg/s]
  REAL(wp)                              :: grav_salt         !<  salt flux moved by gravity drainage between two outputs [kg*ppt/s]
  REAL(wp)                              :: grav_temp         !<  average temperature of gravity drainage brine between two outputs [T]

  !----Flushing
  REAL(wp)                              :: melt_thick        !<  thickness of fully liquid part of top layer [m] 


  !----Lab fluxes
  REAL(wp)                              :: alpha_flux_instable      !<  Proportionality constant which determines energy flux by the temperature difference T_top>T2m [W/C]
  REAL(wp)                              :: alpha_flux_stable        !<  Proportionality constant which determines energy flux by the temperature difference T_top<T2m [W/C]

  !----Flags
  INTEGER :: atmoflux_flag     !< 1: Use mean climatology of Notz, 2: Use imported reanalysis data, 3: use fixed values defined in mo_init 
  INTEGER :: grav_flag         !< 1: no gravity drainage,  2: Gravity drainage, 3: Simple Drainage
  INTEGER :: prescribe_flag    !< 1: nothing happens, 2: prescribed Salinity profile is prescribed at each timestep (does not disable brine dynamics, just overwrites the salinity!)
  INTEGER :: grav_heat_flag    !< 1: nothing happens, 2: compensates heatfluxes in grav_flag = 2
  INTEGER :: flush_heat_flag   !< 1: nothing happens, 2: compensates heatfluxes in flush_flag = 5
  INTEGER :: turb_flag         !< 1: No bottom turbulence, 2: Bottom mixing
  INTEGER :: salt_flag         !< 1: Sea salt, 2: NaCL 
  INTEGER :: boundflux_flag    !< 1: top and bottom cooling plate, 2:top Notz fluxes, bottom cooling plate 3: top flux=a*(T-T_s)
  INTEGER :: flush_flag        !< 1: no flushing, 4:meltwater is removed artificially, 5:vert and horiz flushing, 6: simplified
  INTEGER :: flood_flag        !< 1: no flooding, 2:normal flooding, 3:simple flooding
  INTEGER :: bottom_flag       !< 1: nothing changes, 2: deactivates all bottom layer dynamics, useful for some debugging and idealized tests
  INTEGER :: debug_flag        !< 1: no raw layer output, 2: each layer  is output at every timestep (warning, file size can be very large)
  INTEGER :: precip_flag       !< 0: solid and liquid precipitation, 1:phase determined by T2m 
  INTEGER :: harmonic_flag     !< 1: minimal permeability is used to calculate Rayleigh number, 2:harmonic mean is used for Rayleigh number 
  INTEGER :: tank_flag         !< 1: nothing, 2: S_bu_bottom and bgc_bottom are calculated as if the experiment is conducted in a tank
  INTEGER :: albedo_flag       !< 1: simple albedo, 2: normal albedo, see func_albedo for details 

  
  !##########################################################################################
  !Variables used to import data
  !##########################################################################################
  INTEGER                               :: Length_Input  !< Sets the input length for atmoflux_flag==2, common value of 13169
  REAL(wp), DIMENSION(8280)             :: Tinput        !< used to read in top temperature for field experiment tests, dimension needs to be set in the code
  REAL(wp), DIMENSION(:), ALLOCATABLE   :: fl_sw_input   !< Used to read in sw fluxes from ERA for atmoflux_flag==2
  REAL(wp), DIMENSION(:), ALLOCATABLE   :: fl_lw_input   !< Used to read in lw fluxes from ERA for atmoflux_flag==2
  REAL(wp), DIMENSION(:), ALLOCATABLE   :: T2m_input     !< Used to read in 2Tm from ERA       for atmoflux_flag==2
  REAL(wp), DIMENSION(:), ALLOCATABLE   :: precip_input  !< Used to read in precipitation from ERA for atmoflux_flag==2
  REAL(wp), DIMENSION(:), ALLOCATABLE   :: time_input    !< Used to read in time from ERA for atmoflux_flag==2
  INTEGER                               :: time_counter  !< Keeps track of input data


  !##########################################################################################
  !Chemicals Baby!
  !All arrays needed to support bigogeochemical tracers
  !Chemical matrixes: index 1 defines the chemical, index 2 the layer
  !##########################################################################################
  INTEGER                               :: bgc_flag      !< 1: no bgc, 2:bgc
  INTEGER                               :: N_bgc         !< Number of chemicals
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: fl_brine_bgc  !< Brine fluxes in a matrix, [kg/s], first index is the layer of origin,  and the second index is the layer of arrival
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: bgc_abs       !< Absolute amount of chemicals [kmol] for each tracer
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: bgc_bu        !< Bulk amounts of chemicals [kmol/kg]
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: bgc_br        !< Brine concentrations of chems [kmol/kg]
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: bgc_bottom    !< Bulk concentrations of chems below the ice [kmol/kg]


  !##########################################################################################
  !Variables needed for tank experiments in which concentrations below the ice change over time
  !##########################################################################################
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: bgc_total     !< Total of chems, for lab experiments with a fixed total amount
  REAL(wp)                              :: m_total       !< Total initial water mass, for lab experiments with a fixed total amount
  REAL(wp)                              :: S_total       !< Total initial salt mass, for lab experiments with a fixed total amount
  REAL(wp)                              :: tank_depth    !< water depth in meters, used to calculate concentrations below ice for tank experiments


  
END MODULE mo_data

