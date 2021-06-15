// In principle, this program works just like "scatterplot_selector.cc", but it
// allows A LOT if customisation. First, a customized set of energies or momenta
// can be defined ( meaning that it's not limited anymore to { 40, 60, 120 } eV or
// { 7.5, 8.0, 8.5 } keV/c ). The relevant quantity - energy or momentum - can be
// chosen via the ONLY_ENERGY_SELECTION switch. It can then select if to use
// forward hemisphere only, or both, via the ONLY_FORWARD_HEMISPHERE switch.
// Additional cuts have been implemented via the RADIAL_CUT and TOF_CUT, 
// which apply cuts on splat_radius and TOF. These cuts have been selected so that
// the plong and prad can be successfully reconstructed ( that is, we are far from
// nodes, and from places where the plong/prad vs tof plots overlap ).
// In addition to this, when run, the program will ask the user to input the two
// quantities she wants to plot: the defaults are SR vs TOF, but the quantity can
// be chosen in real time, without the need to modify the source code ( input from
// the user is always needed, even if just to confirm the default quantities ).
// Finally, the greatest improvement over "scatterplot_selector.cc" is that it can
// run the analysis on more than one batch ( up to nine ) of runs, each one chosen
// via their univocal run_id selector ( present also in the name of the .root files
// containing the data ), that identifies runs made with a common configuration
// ( each number is generated when converting the SIMION data to .root files ). The
// runs to be analyzed are to be specified in the run_list array.

{ 
  gROOT->ProcessLine( "gErrorIgnoreLevel = 3001;" );
  gStyle->SetOptStat( 0 ); // Draw histogram without statbox
  gStyle->SetPalette( kRainBow );
  
  // DO NOT COMMENT
  bool ON_TARGET                = false;
  bool RADIAL_CUT               = false;
  bool TOF_CUT                  = false;
  bool FULL_NODES               = false;
  bool VALID_PXYZ               = false;
  bool FORW_HEMI                = false;
  bool BACK_HEMI                = false;
  bool PRODUCE_FILE             = false;
  
  // Run list array
//   int run_list[] = { 115, 215, 315, 415, 515, 615 }; // Standard 8 Gauss BField
//   int run_list[] = { 130, 230, 330, 430, 530, 630 }; // Standard 8 Gauss BField
//   int run_list[] = { 160, 260, 360, 460, 560, 660 }; // Standard 8 Gauss BField
  int run_list[] = { 0 }; // Standard 8 Gauss BField
  
  std::string cutselector;
  cout << endl;
  cout << "Enter a 1 or a 0 if you want to activate a specific cut (default is always 1) [es. 10011]" << endl;
  bool valid_cuts = false; // bool to finally exit the cuts selection
  while( valid_cuts == false ) { 
    valid_cuts = true; // temporarily true for the final productory
    cutselector = "";
    bool valid_los[7];
    while( cutselector.length() != sizeof( valid_los )/sizeof( bool ) ){ 
      cout << "\e[0;33mAugers\e[0m on target : \e[0;33mAugers\e[0m radial cut : \e[0;33mAugers\e[0m TOF cut : \e[0;33mAugers\e[0m full nodes cut : \e[0;33mAugers\e[0m valid Pxyz : \e[0;33mAugers\e[0m in forward hemisphere : \e[0;33mAugers\e[0m in backward hemisphere - ";
      std::getline( std::cin, cutselector );
      if ( cutselector.empty() ) { cutselector = "1111111"; } // used to be 0000000
    }
    if( ( int )cutselector[0] -'0' == 1 ) { ON_TARGET   = true; valid_los[0] = true; } else if ( ( int )cutselector[0] -'0' == 0 ) { ON_TARGET   = false; valid_los[0] = true; } else { valid_los[0] = false; } // only 0's or 1's accepted
    if( ( int )cutselector[1] -'0' == 1 ) { RADIAL_CUT  = true; valid_los[1] = true; } else if ( ( int )cutselector[1] -'0' == 0 ) { RADIAL_CUT  = false; valid_los[1] = true; } else { valid_los[1] = false; } // only 0's or 1's accepted
    if( ( int )cutselector[2] -'0' == 1 ) { TOF_CUT     = true; valid_los[2] = true; } else if ( ( int )cutselector[2] -'0' == 0 ) { TOF_CUT     = false; valid_los[2] = true; } else { valid_los[2] = false; } // only 0's or 1's accepted
    if( ( int )cutselector[3] -'0' == 1 ) { FULL_NODES  = true; valid_los[3] = true; } else if ( ( int )cutselector[3] -'0' == 0 ) { FULL_NODES  = false; valid_los[3] = true; } else { valid_los[3] = false; } // only 0's or 1's accepted
    if( ( int )cutselector[4] -'0' == 1 ) { VALID_PXYZ  = true; valid_los[4] = true; } else if ( ( int )cutselector[4] -'0' == 0 ) { VALID_PXYZ  = false; valid_los[4] = true; } else { valid_los[4] = false; } // only 0's or 1's accepted
    if( ( int )cutselector[5] -'0' == 1 ) { FORW_HEMI   = true; valid_los[5] = true; } else if ( ( int )cutselector[5] -'0' == 0 ) { FORW_HEMI   = false; valid_los[5] = true; } else { valid_los[5] = false; } // only 0's or 1's accepted
    if( ( int )cutselector[6] -'0' == 1 ) { BACK_HEMI   = true; valid_los[6] = true; } else if ( ( int )cutselector[6] -'0' == 0 ) { BACK_HEMI   = false; valid_los[6] = true; } else { valid_los[6] = false; } // only 0's or 1's accepted
    for( int i=0;i<sizeof( valid_los )/sizeof( bool );++i ){ valid_cuts *= valid_los[i]; } // Check on the validity of the entered string
  }
  
  // Revamped TString structure
  // This is a matrix of TStrings containing 3 columns and a number of rows equal to the number of analyses possible
  // The first column contains a plain text definition of a certain quantity;
  // the second column contains a LaTeX string, for the purposes of naming histograms;
  // the third column contains the associated unit of measurement
  // the fourth column contains data to title an eventual .root output file
  TString outputs[][4] = { 
    { "TOF"                                                                         , "TOF"                                                           , "(#mus)"   , "TOF"                                             }, //   0
    { "splat_radius"                                                                , "splat_radius"                                                  , "(mm)"     , "splat_radius"                                    }, //   1
    { "Simulated elevation"                                                         , "init_elevation"                                                , "(deg)"    , "init_elevation"                                  }, //   2
    { "splat_azimuth"                                                               , "splat_azimuth"                                                 , "(deg)"    , "splat_azimuth"                                   }, //   3
    { "init kinetic energy"                                                         , "init kinetic energy"                                           , "(eV)"     , "init_ke"                                         }, //   4
    { "Simulated prad"                                                              , "init_p_{rad}"                                                  , "(keV/c)"  , "init_prad"                                       }, //   5
    { "Simulated plong"                                                             , "init_p_{long}"                                                 , "(keV/c)"  , "init_plong"                                      }, //   6
    { "Simulated ptot"                                                              , "init_p_{tot}"                                                  , "(keV/c)"  , "init_ptot"                                       }, //   7
    { "Simulated azimuth"                                                           , "init_#varphi"                                                  , "(deg)"    , "init_azimuth"                                    }, //   8
    { "Simulated px"                                                                , "init_p_{x}"                                                    , "(keV/c)"  , "init_px"                                         }, //   9
    { "Simulated py"                                                                , "init_p_{y}"                                                    , "(keV/c)"  , "init_py"                                         }, //  10
    { "Simulated pz"                                                                , "init_p_{z}"                                                    , "(keV/c)"  , "init_pz"                                         }, //  11
    { "(iso) Uncertainty on calculated kinetic energy"                              , "#delta (iso) calc kinetic energy"                              , "(eV)"     , "delta_calc_ke"                                   }, //  12
    { "(iso) Uncertainty on calculated p_rad"                                       , "#delta (iso) calc p_{rad}"                                     , "(keV/c)"  , "delta_calc_prad"                                 }, //  13
    { "(iso) Uncertainty on calculated p_long"                                      , "#delta (iso) calc p_{long}"                                    , "(keV/c)"  , "delta_calc_plong"                                }, //  14
    { "(iso) Uncertainty on calculated p_tot"                                       , "#delta (iso) calc p_{tot}"                                     , "(keV/c)"  , "delta_calc_ptot"                                 }, //  15
    { "(iso) Uncertainty on calculated azimuth"                                     , "#delta (iso) calc #varphi"                                     , "(deg)"    , "delta_calc_azimuth"                              }, //  16
    { "(iso) Uncertainty on calculated elevation"                                   , "#delta (iso) calc #theta"                                      , "(deg)"    , "delta_calc_elevation"                            }, //  17
    { "(iso) Calculated kinetic energy"                                             , "(iso) calc kinetic energy"                                     , "(eV)"     , "calc_ke"                                         }, //  18
    { "(iso) Calculated p_rad"                                                      , "(iso) calc p_{rad}"                                            , "(keV/c)"  , "calc_prad"                                       }, //  19
    { "(iso) Calculated p_long"                                                     , "(iso) calc p_{long}"                                           , "(keV/c)"  , "calc_plong"                                      }, //  20
    { "(iso) Calculated p_tot"                                                      , "(iso) calc p_{tot}"                                            , "(keV/c)"  , "calc_ptot"                                       }, //  21
    { "(iso) Calculated azimuth"                                                    , "(iso) calc #varphi"                                            , "(deg)"    , "calc_azimuth"                                    }, //  22
    { "(iso) Calculated elevation"                                                  , "(iso) calc #theta"                                             , "(deg)"    , "calc_elevation"                                  }, //  23
    { "(iso) Calculated p_x"                                                        , "(iso) calc p_{x}"                                              , "(keV/c)"  , "calc_px"                                         }, //  24
    { "(iso) Calculated p_y"                                                        , "(iso) calc p_{y}"                                              , "(keV/c)"  , "calc_py"                                         }, //  25
    { "(iso) Calculated p_z"                                                        , "(iso) calc p_{z}"                                              , "(keV/c)"  , "calc_pz"                                         }, //  26
    { "(COLTRIMS) Theoretical p_rad"                                                , "(COLTRIMS) theo p_{rad}"                                       , "(keV/c)"  , "theo_prad"                                       }, //  27
    { "(COLTRIMS) Theoretical p_long"                                               , "(COLTRIMS) theo p_{long}"                                      , "(keV/c)"  , "theo_plong"                                      }, //  28
    { "(COLTRIMS) Theoretical azimuth"                                              , "(COLTRIMS) theo #varphi"                                       , "(deg)"    , "theo_azimuth"                                    }, //  29
    { "(COLTRIMS) Theoretical elevation"                                            , "(COLTRIMS) theo #theta"                                        , "(deg)"    , "theo_elevation"                                  }, //  30
    { "Uncertainty on (COLTRIMS) theoretical p_rad"                                 , "#delta (COLTRIMS) theo p_{rad}"                                , "(keV/c)"  , "delta_theo_prad"                                 }, //  31
    { "Uncertainty on (COLTRIMS) theoretical p_long"                                , "#delta (COLTRIMS) theo p_{long}"                               , "(keV/c)"  , "delta_theo_plong"                                }, //  32
    { "Uncertainty on (COLTRIMS) theoretical p_tot"                                 , "#delta (COLTRIMS) theo p_{tot}"                                , "(keV/c)"  , "delta_theo_ptot"                                 }, //  33
    { "Uncertainty on (COLTRIMS) theoretical azimuth"                               , "#delta (COLTRIMS) theo #varphi"                                , "(deg)"    , "delta_theo_azimuth"                              }, //  34
    { "Uncertainty on (COLTRIMS) theoretical elevation"                             , "#delta (COLTRIMS) theo #theta"                                 , "(deg)"    , "delta_theo_elevation"                            }, //  35
    { "[ Uncertainty on (iso) calculated p_rad ] / simulated p_rad"                 , "#delta (iso) calc p_{rad} / init p_{rad}"                      , ""         , "delta_calc_prad_over_init_prad"                  }, //  36
    { "[ Uncertainty on (iso) calculated p_long ] / simulated p_long"               , "#delta (iso) calc p_{long} / init p_{long}"                    , ""         , "delta_calc_plong_over_init_plong"                }, //  37
    { "[ Uncertainty on (iso) calculated p_tot ] / simulated p_tot"                 , "#delta (iso) calc p_{tot} / init p_{tot}"                      , ""         , "delta_calc_ptot_over_init_ptot"                  }, //  38
    { "[ Uncertainty on (iso) calculated p_rad ] / (iso) calculated p_rad"          , "#delta (iso) calc p_{rad} / (iso) calc p_{rad}"                , ""         , "delta_calc_prad_over_calc_prad"                  }, //  39
    { "[ Uncertainty on (iso) calculated p_long ] / (iso) calculated p_long"        , "#delta (iso) calc p_{long} / (iso) calc p_{long}"              , ""         , "delta_calc_plong_over_calc_plong"                }, //  40
    { "[ Uncertainty on (iso) calculated p_tot ] / (iso) calculated p_tot"          , "#delta (iso) calc p_{tot} / (iso) calc p_{tot}"                , ""         , "delta_calc_ptot_over_calc_ptot"                  }, //  41
    { "(iso) calculated p_rad - simulated p_rad"                                    , "(iso) calc p_{rad} - init p_{rad}"                             , "(keV/c)"  , "calc_prad_minus_init_prad"                       }, //  42
    
    { "[ (iso) calculated p_rad - simulated p_rad ] / simulated p_rad"              , "[ (iso) calc p_{rad} - init p_{rad} ] / init p_{rad}"          , ""         , "calc_prad_minus_init_prad_over_init_prad"        }, //  43
    { "(iso) calculated p_long - simulated p_long"                                  , "(iso) calc p_{long} - init p_{long}"                           , "(keV/c)"  , "calc_plong_minus_init_plong"                     }, //  44
    { "[ (iso) calculated p_long - simulated p_long ] / simulated p_long"           , "[ (iso) calc p_{long} - init p_{long} ] / init p_{long}"       , ""         , "calc_plong_minus_init_plong_over_init_plong"     }, //  45
    { "(iso) calculated p_tot - simulated p_tot"                                    , "(iso) calc p_{tot} - init p_{tot}"                             , "(keV/c)"  , "calc_ptot_minus_init_ptot"                       }, //  46
    { "[ (iso) calculated p_tot - simulated p_tot ] / simulated p_tot"              , "[ (iso) calc p_{tot} - init p_{tot} ] / init p_{tot}"          , ""         , "calc_ptot_minus_init_ptot_over_init_ptot"        }, //  47
    { "(iso) calculated azimuth - simulated azimuth"                                , "(iso) calc #varphi - init #varphi"                             , "(deg)"    , "calc_azimuth_minus_init_azimuth"                 }, //  48
    { "(iso) calculated elevation - simulated elevation"                            , "(iso) calc #theta - init #theta"                               , "(deg)"    , "calc_elevation_minus_init_elevation"             }, //  49
    { "(iso) calculated kinetic energy - simulated kinetic energy"                  , "(iso) calc kinetic energy - init kinetic energy"               , "(eV)"     , "calc_ke_minus_init_ke"                           }, //  50
    { "(COLTRIMS) theoretical p_rad - simulated p_rad"                              , "(COLTRIMS) theo p_{rad} - init p_{rad}"                        , "(keV/c)"  , "theo_prad_minus_init_prad"                       }, //  51
    { "[ (COLTRIMS) theoretical p_rad - simulated p_rad ] / simulated p_rad"        , "[ (COLTRIMS) theo p_{rad} - init p_{rad} ] / init p_{rad}"     , ""         , "theo_prad_minus_init_prad_over_init_prad"        }, //  52
    { "(COLTRIMS) theoretical p_long - simulated p_long"                            , "(COLTRIMS) theo p_{long} - init p_{long}"                      , "(keV/c)"  , "theo_plong_minus_init_plong"                     }, //  53
    { "[ (COLTRIMS) theoretical p_long - simulated p_long ] / simulated p_long"     , "[ (COLTRIMS) theo p_{long} - init p_{long} ] / init p_{long}"  , ""         , "theo_plong_minus_init_plong_over_init_plong"     }, //  54
    { "(COLTRIMS) theoretical p_tot - simulated p_tot"                              , "(COLTRIMS) theo p_{tot} - init p_{tot}"                        , "(keV/c)"  , "theo_ptot_minus_init_ptot"                       }, //  55
    { "[ (COLTRIMS) theoretical p_tot - simulated p_tot ] / simulated p_tot"        , "[ (COLTRIMS) theo p_{tot} - init p_{tot} ] / init p_{tot}"     , ""         , "theo_ptot_minus_init_ptot_over_init_ptot"        }, //  56
    { "(COLTRIMS) theoretical azimuth - simulated azimuth"                          , "(COLTRIMS) theo #varphi - init #varphi"                        , "(deg)"    , "theo_azimuth_minus_azimuth"                      }, //  57
    { "(COLTRIMS) theoretical elevation - simulated elevation"                      , "(COLTRIMS) theo #theta - init #theta"                          , "(deg)"    , "theo_elevation_minus_elevation"                  }, //  58
    { "(COLTRIMS) calculated kinetic energy - simulated kinetic energy"             , "(COLTRIMS) calc kinetic energy - init kinetic energy"          , "(eV)"     , "theo_ke_minus_init_ke"                           }, //  59
    { "( wt - ⌊( wt/2π )⌋·2π )/2"                                                   , "( #omegat - floor[( #omegat/2#pi )]#upoint2#pi )/2"            , ""         , "wt_floor"                                        }, //  60
    { "sin( wt/2 )"                                                                 , "sin( #omega#upointt/2 )"                                       , ""         , "sin_wt2"                                         }, //  61
    { "|sin( wt/2 )|"                                                               , "|sin( #omega#upointt/2 )|"                                     , ""         , "abs_sin_wt2"                                     }, //  62
    { "splat_x"                                                                     , "splat_x"                                                       , "(mm)"     , "splat_x"                                         }, //  63
    { "splat_y"                                                                     , "splat_y"                                                       , "(mm)"     , "splat_y"                                         }, //  64
    { "splat_z"                                                                     , "splat_z"                                                       , "(mm)"     , "splat_z"                                         }, //  65
    { "Simulated x"                                                                 , "Simulated x"                                                   , "(mm)"     , "init_x"                                          }, //  66
    { "Simulated y"                                                                 , "Simulated y"                                                   , "(mm)"     , "init_y"                                          }, //  67
    { "Simulated z"                                                                 , "Simulated z"                                                   , "(mm)"     , "init_z"                                          }, //  68
    { "(iso) calculated p_rad - (COLTRIMS) theoretical p_rad"                       , "(iso) calc p_{rad} - (COLTRIMS) theo p_{rad}"                  , "(keV/c)"  , "calc_prad_minus_theo_prad"                       }, //  69
    { "(iso) calculated p_long - (COLTRIMS) theoretical p_long"                     , "(iso) calc p_{long} - (COLTRIMS) theo p_{long}"                , "(keV/c)"  , "calc_plong_minus_theo_plong"                     }, //  70
    { "(iso) calculated p_tot - (COLTRIMS) theoretical p_tot"                       , "(iso) calc p_{tot} - (COLTRIMS) theo p_{tot}"                  , "(keV/c)"  , "calc_ptot_minus_theo_ptot"                       }, //  71
    { "(iso) calculated azimuth - (COLTRIMS) theoretical azimuth"                   , "(iso) calc azimuth - (COLTRIMS) theo azimuth"                  , "(deg)"    , "calc_azimuth_minus_theo_azimuth"                 }, //  72
    { "(iso) calculated elevation - (COLTRIMS) theoretical elevation"               , "(iso) calc elevation - (COLTRIMS) theo elevation"              , "(deg)"    , "calc_elevation_minus_theo_elevation"             }, //  73
    { "(iso) calculated kinetic energy - (COLTRIMS) theoretical kinetic energy"     , "(iso) calc kinetic energy - (COLTRIMS) theo kinetic energy"    , "(eV)"     , "calc_ke_minus_theo_ke"                           }, //  74 MCP MEASURES FROM HERE
    { "[MCP] Uncertainty on (iso) calculated kinetic energy"                        , "[MCP] #delta (iso) calc kinetic energy"                        , "(eV)"     , "mcp_delta_calc_ke"                               }, //  75
    { "[MCP] Uncertainty on (iso) calculated p_rad"                                 , "[MCP] #delta (iso) calc p_{rad}"                               , "(keV/c)"  , "mcp_delta_calc_prad"                             }, //  76
    { "[MCP] Uncertainty on (iso) calculated p_long"                                , "[MCP] #delta (iso) calc p_{long}"                              , "(keV/c)"  , "mcp_delta_calc_plong"                            }, //  77
    { "[MCP] Uncertainty on (iso) calculated p_tot"                                 , "[MCP] #delta (iso) calc p_{tot}"                               , "(keV/c)"  , "mcp_delta_calc_ptot"                             }, //  78
    { "[MCP] Uncertainty on (iso) calculated azimuth"                               , "[MCP] #delta (iso) calc #varphi"                               , "(deg)"    , "mcp_delta_calc_azimuth"                          }, //  79
    { "[MCP] Uncertainty on (iso) calculated elevation"                             , "[MCP] #delta (iso) calc #theta"                                , "(deg)"    , "mcp_delta_calc_elevation"                        }, //  80
    { "[MCP] (iso) Calculated kinetic energy"                                       , "[MCP] (iso) calc kinetic energy"                               , "(eV)"     , "mcp_calc_ke"                                     }, //  81
    { "[MCP] (iso) Calculated p_rad"                                                , "[MCP] (iso) calc p_{rad}"                                      , "(keV/c)"  , "mcp_calc_prad"                                   }, //  82
    { "[MCP] (iso) Calculated p_long"                                               , "[MCP] (iso) calc p_{long}"                                     , "(keV/c)"  , "mcp_calc_plong"                                  }, //  83
    { "[MCP] (iso) Calculated p_tot"                                                , "[MCP] (iso) calc p_{tot}"                                      , "(keV/c)"  , "mcp_calc_ptot"                                   }, //  84
    { "[MCP] (iso) Calculated azimuth"                                              , "[MCP] (iso) calc #varphi"                                      , "(deg)"    , "mcp_calc_azimuth"                                }, //  85
    { "[MCP] (iso) Calculated elevation"                                            , "[MCP] (iso) calc #theta"                                       , "(deg)"    , "mcp_calc_elevation"                              }, //  86
    { "[MCP] (iso) Calculated p_x"                                                  , "[MCP] (iso) calc p_{x}"                                        , "(keV/c)"  , "mcp_calc_px"                                     }, //  87
    { "[MCP] (iso) Calculated p_y"                                                  , "[MCP] (iso) calc p_{y}"                                        , "(keV/c)"  , "mcp_calc_py"                                     }, //  88
    { "[MCP] (iso) Calculated p_z"                                                  , "[MCP] (iso) calc p_{z}"                                        , "(keV/c)"  , "mcp_calc_pz"                                     }, //  89
    { "[MCP] [ Uncertainty on (iso) calculated p_rad ] / simulated p_rad"           , "[MCP] #delta (iso) calc p_{rad} / init p_{rad}"                , ""         , "mcp_delta_calc_prad_over_init_prad"              }, //  90
    { "[MCP] [ Uncertainty on (iso) calculated p_long ] / simulated p_long"         , "[MCP] #delta (iso) calc p_{long} / init p_{long}"              , ""         , "mcp_delta_calc_plong_over_init_plong"            }, //  91
    { "[MCP] [ Uncertainty on (iso) calculated p_tot ] / simulated p_tot"           , "[MCP] #delta (iso) calc p_{tot} / init p_{tot}"                , ""         , "mcp_delta_calc_ptot_over_init_ptot"              }, //  92
    { "[MCP] [ Uncertainty on (iso) calculated p_rad ] / (iso) calculated p_rad"    , "[MCP] #delta (iso) calc p_{rad} / (iso) calc p_{rad}"          , ""         , "mcp_delta_calc_prad_over_calc_prad"              }, //  93
    { "[MCP] [ Uncertainty on (iso) calculated p_long ] / (iso) calculated p_long"  , "[MCP] #delta (iso) calc p_{long} / (iso) calc p_{long}"        , ""         , "mcp_delta_calc_plong_over_calc_plong"            }, //  94
    { "[MCP] [ Uncertainty on (iso) calculated p_tot ] / (iso) calculated p_tot"    , "[MCP] #delta (iso) calc p_{tot} / (iso) calc p_{tot}"          , ""         , "mcp_delta_calc_ptot_over_calc_ptot"              }, //  95
    { "[MCP] (iso) calculated p_rad - simulated p_rad"                              , "[MCP] (iso) calc p_{rad} - init p_{rad}"                       , "(keV/c)"  , "mcp_calc_prad_minus_init_prad"                   }, //  96
    { "[MCP] [ (iso) calculated p_rad - simulated p_rad ] / simulated p_rad"        , "[MCP] [ (iso) calc p_{rad} - init p_{rad} ] / init p_{rad}"    , ""         , "mcp_calc_prad_minus_init_prad_over_init_prad"    }, //  97
    { "[MCP] (iso) calculated p_long - simulated p_long"                            , "[MCP] (iso) calc p_{long} - init p_{long}"                     , "(keV/c)"  , "mcp_calc_plong_minus_init_plong"                 }, //  98
    { "[MCP] [ (iso) calculated p_long - simulated p_long ] / simulated p_long"     , "[MCP] [ (iso) calc p_{long} - init p_{long} ] / init p_{long}" , ""         , "mcp_calc_plong_minus_init_plong_over_init_plong" }, //  99
    { "[MCP] (iso) calculated p_tot - simulated p_tot"                              , "[MCP] (iso) calc p_{tot} - init p_{tot}"                       , "(keV/c)"  , "mcp_calc_ptot_minus_init_ptot"                   }, // 100
    { "[MCP] [ (iso) calculated p_tot - simulated p_tot ] / simulated p_tot"        , "[MCP] [ (iso) calc p_{tot} - init p_{tot} ] / init p_{tot}"    , ""         , "mcp_calc_ptot_minus_init_ptot_over_init_ptot"    }, // 101
    { "[MCP] (iso) calculated azimuth - simulated azimuth"                          , "[MCP] (iso) calc #varphi - init #varphi"                       , "(deg)"    , "mcp_calc_azimuth_minus_init_azimuth"             }, // 102
    { "[MCP] (iso) calculated elevation - simulated elevation"                      , "[MCP] (iso) calc #theta - init #theta"                         , "(deg)"    , "mcp_calc_elevation_minus_init_elevation"         }, // 103
    { "[MCP] (iso) calculated kinetic energy - simulated kinetic energy"            , "[MCP] (iso) calc kinetic energy - init kinetic energy"         , "(eV)"     , "mcp_calc_ke_minus_init_ke"                       },  // 104
    { "(iso) interpolated p_rad - simulated p_rad"                                  , "(iso) intr p_{rad} - init p_{rad}"                             , "(keV/c)"  , "inte_prad_minus_init_prad"                       } //  105
  };
  
  // Number of total analyses ( must be divided by 3 because outputs is a Nx3 matrix )
  int nanalyses = sizeof( outputs )/sizeof( TString )/4;
  
  // List of available analyses
  cout << endl;
  cout << "ANALYSIS SELECTOR" << endl;
  for( int i=0;i<nanalyses;++i ){ cout << "  case " << i << "\t: " << outputs[i][0].Data() << endl; }
  cout << endl;
  
  int select_x = 0;
  int select_y = 1;
  int select_z = 4;
  
  // Vector that contains a pair<selection, axis>
  vector< pair<int, int> > *vpin = new vector< pair<int, int> >;
  
  // Dummy string for the selection of the quantities to plot
  std::string input;
  
  // Selection of x quantity
  cout << "Select x axis (default '0 : TOF'): ";
  std::getline( std::cin, input );
  if ( !input.empty() ) { 
    std::istringstream stream( input );
    stream >> select_x;
  }
  vpin->push_back( make_pair( select_x, 0 ) );
  cout << outputs[select_x][0] << " selected" << endl;
  cout << endl;
  
  // Selection of y quantity
  cout << "Select y axis (default '1 : splat_radius'): ";
  std::getline( std::cin, input );
  if ( !input.empty() ) { 
    std::istringstream stream( input );
    stream >> select_y;
  }
  vpin->push_back( make_pair( select_y, 1 ) );
  cout << outputs[select_y][0] << " selected" << endl;
  cout << endl;
  
  // Selection of z quantity
  cout << "Select z axis (default '4 : init kinetic energy'): ";
  std::getline( std::cin, input );
  if ( !input.empty() ) { 
    std::istringstream stream( input );
    stream >> select_z;
  }
  vpin->push_back( make_pair( select_z, 2 ) );
  cout << outputs[select_z][0] << " selected" << endl;
  cout << endl;
  
  // Sorting of vpin ( explaination follows )
  sort( vpin->begin(), vpin->end() );
  // EXPLAINATION: vpin contains three pairs<selection, axis>, and it's sorted about
  // the selection. The program will then initialize a three-element array which will
  // contain, at the i-th position, the value matched to the SECOND element of the
  // pair ( the one associated to the axis ).
  //
  // EXAMPLE:
  //   vpin = ( <28, 0> ; <15, 1> ; <40, 2> ).
  // After sorting it's:
  //   vpin = ( <15, 1> ; <28, 0> ; <40, 2> ).
  // There will be an array vpout[3], which will be filled via:
  //   if( vpin->at(z).first ==  NN ){ vpout[vpin->at(z).second] = 'the quantity associated to NN' }.
  // For this example, vpout will be filled as:
  //   vpin->at( 0 ).first == 15  =>  vpout[vpin->at( 0 ).second] = vpout[1] = 'the quantity associated to 15';
  //   vpin->at( 1 ).first == 28  =>  vpout[vpin->at( 1 ).second] = vpout[0] = 'the quantity associated to 28';
  //   vpin->at( 2 ).first == 40  =>  vpout[vpin->at( 2 ).second] = vpout[2] = 'the quantity associated to 40';
  // to give:
  //   vpout[0] = 'the quantity associated to 28';
  //   vpout[1] = 'the quantity associated to 15';
  //   vpout[2] = 'the quantity associated to 40';
  // which is the initial order!!!
  
  int nruns = sizeof( run_list ) / sizeof( int );
  
  TString cuts = "init_ke > 0.";
  //   cuts += "&& ( ( splat_radius >=  9.95 && splat_radius <= 10.05 )";
  //   cuts += "|| ( splat_radius >= 19.95 && splat_radius <= 20.05 )";
  //   cuts += "|| ( splat_radius >= 29.95 && splat_radius <= 30.05 )";
  //   cuts += "|| ( splat_radius >= 39.95 && splat_radius <= 40.05 )";
  //   cuts += "|| ( splat_radius >= 49.95 && splat_radius <= 50.05 ) )";
  if( ON_TARGET )       { cuts += TString::Format( " && splat_radius >= 0 && splat_radius <= %.2f", CX_SR_HIGH ); }
  if( RADIAL_CUT ) { cuts += TString::Format( " && splat_radius >= %.2f && splat_radius <= %.2f", CX_SR_LOW, CX_SR_HIGH ); }
  //   if( TOF_CUT )    { cuts += TString::Format( " && ( ( init_ke > 90. && ( ( tof < %.8f ) || ( tof > %.8f && tof < %.8f ) ) ) || ( init_ke < 90. && tof <= 0.18 ) || ( init_ke < 90. && tof > 0.18 && tof < %.8f && splat_radius > ( %f + %f*tof + %f*tof*tof + %f*tof*tof*tof ) ) )", CX_TOF_HI, CX_TOF_HI_LOW, CX_TOF_HI_HIGH, CX_TOF_LO, CX_TOF_p0, CX_TOF_p1, CX_TOF_p2, CX_TOF_p3 ); }
  if( TOF_CUT )    { cuts += TString::Format( " && tof < %.8f ", CX_TOF ); }
  //   if( TOF_CUT )    { cuts += TString::Format( " && ( ( init_ke > 90. && tof < %.8f ) || ( init_ke < 90. && tof <= 0.18 ) )", CX_TOF_HI ); }
  //   if( TOF_CUT )    { cuts += TString::Format( " && tof < %.4f", CX_TOF ); }
  if( FULL_NODES ) { 
    for( int nnodes=0; nnodes<sizeof( nodes )/sizeof( double ); ++nnodes ){ cuts += TString::Format( " && ( tof < %.4f || tof > %.4f )", nodes[nnodes]-CX_NODE, nodes[nnodes]+CX_NODE ); }
  }
  //   if( VALID_PXYZ )       { cuts += " && calc_prad() != 0. && calc_plong() != 0."; }  // NON FUNZIONA!!!
  if( FORW_HEMI ^ BACK_HEMI ){ 
    if( FORW_HEMI )      { cuts += " && TMath::ACos( init_dircos )*TMath::RadToDeg() >= 0. && TMath::ACos( init_dircos )*TMath::RadToDeg() <=  90."; }
    if( BACK_HEMI )      { cuts += " && TMath::ACos( init_dircos )*TMath::RadToDeg() > 90. && TMath::ACos( init_dircos )*TMath::RadToDeg() <= 180."; }
  }
  else                 { cuts += " && TMath::ACos( init_dircos )*TMath::RadToDeg() >= 0. && TMath::ACos( init_dircos )*TMath::RadToDeg() <= 180."; }
  //   cuts += " && init_dircos >= 0.93969 && init_dircos <= 1.";
  
  // Accessory variables
  TString var_x, var_y, var_z;
  TString unit_x, unit_y, unit_z;
  
  // Quantity selectors
  var_x = outputs[select_x][1]; unit_x = outputs[select_x][2];
  var_y = outputs[select_y][1]; unit_y = outputs[select_y][2];
  var_z = outputs[select_z][1]; unit_z = outputs[select_z][2];
  
  // Output file name
  TString outfile_name;
  
  // Selector to produce output .root file
  std::string outfile;
  cout << "Do you want to produce an output .root file [Y/N]? (default 'N') : ";
  std::getline( std::cin, outfile );
  if ( outfile == "y" || outfile == "Y" ) { 
    PRODUCE_FILE = true;
//     outfile_name = TString::Format( "%d-%s-%s-%s-%s", run_list[0], outputs[select_x][3].Data(), outputs[select_y][3].Data(), outputs[select_z][3].Data(), cutselector.c_str() );
    outfile_name = TString::Format( "%s-%s-%s-%s", outputs[select_x][3].Data(), outputs[select_y][3].Data(), outputs[select_z][3].Data(), cutselector.c_str() );    // RUN_ID AGNOSTIC
  }
  cout << endl;
  
  // Titling of the final histogram and axes
  TString title   = TString::Format( "%s vs %s vs %s", var_y.Data(), var_x.Data(), var_z.Data() );
  TString axis_x  = TString::Format( "%s %s", var_x.Data(), unit_x.Data() );
  TString axis_y  = TString::Format( "%s %s", var_y.Data(), unit_y.Data() );
  TString axis_z  = TString::Format( "%s %s", var_z.Data(), unit_z.Data() );
  
  // Recap of the selected analysis
  if( 
    ( ON_TARGET + RADIAL_CUT + TOF_CUT + FULL_NODES + FORW_HEMI + BACK_HEMI )
  )                 { cout << "\e[4;35mUsing additional cuts:\e[0m"                                                                                 << endl; }
  else              { cout << "No cuts applied"                                                                                                     << endl; }
  if( ON_TARGET )   { cout << " - \e[0;33mAugers\e[0m hitting the MCP"                                                                              << endl; }
  if( RADIAL_CUT )  { cout << TString::Format( " - \e[0;33mAugers\e[0m splat_radius cut : %.0f <= splat_radius <= %.0f mm", CX_SR_LOW, CX_SR_HIGH ) << endl; }
  //   if( TOF_CUT )     { cout << TString::Format( " - \e[0;33mAugers\e[0m TOF cut : TOF < %.4f us", CX_TOF )                                            << endl; }
  //   if( TOF_CUT )     { cout << " - \e[0;33mAugers\e[0m over TOF/SR pol2 curves"                                                                     << endl; }
  if( FULL_NODES )  { cout << TString::Format( " - \e[0;33mAugers\e[0m full nodes cut (width of cut = %.4f us)", CX_NODE )                        << endl; }
  if( VALID_PXYZ )  { cout << " - prad and plong successfully reconstructed for \e[0;33mAugers\e[0m"                                                << endl; }
  if( FORW_HEMI )   { cout << " - \e[0;33mAugers\e[0m in forward hemisphere"                                                                        << endl; }
  if( BACK_HEMI )   { cout << " - \e[0;33mAugers\e[0m in backward hemisphere"                                                                       << endl; }
  cout << endl;
  
  cout << "\e[4;35mProducing plot for\e[0m" << endl;
  cout << title.Data() << endl;
  cout << endl;
  
  cout << "\e[4;35mUsing configurations:\e[0m" << endl;
  for( int rid = 0; rid < sizeof( run_list ) / sizeof( int ); ++rid ){ 
    for( int k = 0; k < so->GetListOfFiles()->GetEntries(); ++k ){ 
      TString dummy = ( TString )so->GetListOfFiles()->At( k )->GetTitle();
      if( dummy.SubString( TString::Format( "_%d.", run_list[rid] ) ) == "" ){ continue; }
      cout << dummy.Data() << endl;
    }
  }
  cout << endl;
  
  cout << "\e[4;35mCut configuration\e[0m" << endl;
  cout << cutselector << endl;
  cout << endl;
  
  if( PRODUCE_FILE ){ 
    cout << "\e[4;35mProducing output file:\e[0m \e[93m" << outfile_name.Data() << ".root\e[0m" << endl;
    if( !( gSystem->AccessPathName( TString::Format( "%s.root", outfile_name.Data() ) ) ) ){ cout << "\e[41;4mWARNING: THE FILE ALREADY EXISTS, POSSIBLE OVERWRITING!!!\e[0m" << endl; }
    cout << endl;
  }
  
  // Validating analysis - the user is able to exit if analysis is not correct
  std::string yesno;
  cout << "\e[7mDo you want to continue with this analysis [Y/N]? (default 'Y') :\e[0m ";
  std::getline( std::cin, yesno );
  if ( yesno == "n" || yesno == "N" ) { 
    gApplication->Terminate( 1 );
  }
  cout << endl;
  
  cout << "Number of total events = ";
//   vector<int> *ntracks = FindTrack( TString::Format( "run_id == %d", run_list[0] ) );
  vector<int> *ntracks = FindTrack( TString::Format("") );   // RUN_ID AGNOSTIC
  cout << ntracks->size() << endl;
  
  
  // Output tree - only for working reasons
  TTree *tout = new TTree( "tout", "tout" );
  double outx, outy, outz;
  long long outn;
  tout->Branch( "outx", &outx );
  tout->Branch( "outy", &outy );
  tout->Branch( "outz", &outz );
  tout->Branch( "outn", &outn );
  
  // Time counter
  TTimeStamp *time_counter = new TTimeStamp();
  
  // LONG EXPLAINATION
  // Completely rewritten structure. This only works for a single run_id. It uses the
  // ( very fast ) "Draw" method of TTree to draw 2D scatter plots which have a third
  // axis represented with the color of each point.
  // The vector<int> events are scanned, and each selected quantity is stored into a
  // temporary TTree, which is then used for the "Draw" method.
  
  // Loop on runs
//   for( int run = 0; run < nruns; ++run ){ 
  for( int run = 0; run < 1; ++run ){   // RUN_ID AGNOSTIC
    
    double vpout[3] = { 0., 0., 0. };
    
    // Selection of the required Tracks
//     vector<int> *sng_tracks_tgt = FindTrack( TString::Format( "run_id == %d && %s", run_list[run], cuts.Data() ) );
    vector<int> *sng_tracks_tgt = FindTrack( TString::Format( "%s", cuts.Data() ) );    // RUN_ID AGNOSTIC
    
    cout << "Number of events after cuts = ";
    cout << sng_tracks_tgt->size() << endl;
    cout << endl;
    
    cout << "Filling the histogram" << endl;
    int niter = sng_tracks_tgt->size();
    int step_delta, step;
    if( niter > 10. ){ 
      step_delta = -1; // smaller values ( even < 0 ) = faster counter; larger values = slower counter
      step = TMath::Power( 10, step_delta + TMath::Floor( TMath::Log10( niter ) ) );    // Event counter
    }
    
    time_counter->Set();
    double time_start = time_counter->AsDouble();
    double time_previous = time_start;
    //     int time_start = time_counter->GetSec();
    //     int time_previous = time_start;
    // Filling of of the TGraph2DErrorss
    for( int i=0;i<sng_tracks_tgt->size();++i ){ 
      so->GetEntry( sng_tracks_tgt->at( i ) );
      if( niter > 10. ){ 
        if( ( i )%step == 0 || i == niter-1 ){ 
          time_counter->Set();
          double time_now = time_counter->AsDouble();
          cout << i+1 << "/" << niter << "\t" << TString::Format( "%.2f", time_now - time_start ) << "s elapsed" << "\t Δt = " << TString::Format( "%.2f", time_now - time_previous ) << "s" << "\tEntry #:" << so->GetReadEntry() << endl;
          time_previous = time_now;
        }
      }
      
      for( int z=0; z<3; ++z ){ 
        if( vpin->at(z).first ==   0 ){                                            vpout[vpin->at(z).second] = e.tof;                                                                                                     }
        if( vpin->at(z).first ==   1 ){                                            vpout[vpin->at(z).second] = e.splat_radius;                                                                                            }
        if( vpin->at(z).first ==   2 ){                                            vpout[vpin->at(z).second] = e.init_elevation();                                                                                        }
        if( vpin->at(z).first ==   3 ){                                            vpout[vpin->at(z).second] = e.splat_azimuth();                                                                                         }
        if( vpin->at(z).first ==   4 ){                                            vpout[vpin->at(z).second] = e.init_ke;                                                                                                 }
        if( vpin->at(z).first ==   5 ){                                            vpout[vpin->at(z).second] = e.init_prad;                                                                                               }
        if( vpin->at(z).first ==   6 ){                                            vpout[vpin->at(z).second] = e.init_plong;                                                                                              }
        if( vpin->at(z).first ==   7 ){                                            vpout[vpin->at(z).second] = e.init_ptot;                                                                                               }
        if( vpin->at(z).first ==   8 ){                                            vpout[vpin->at(z).second] = e.init_azimuth();                                                                                          }
        if( vpin->at(z).first ==   9 ){                                            vpout[vpin->at(z).second] = e.init_px;                                                                                                 }
        if( vpin->at(z).first ==  10 ){                                            vpout[vpin->at(z).second] = e.init_py;                                                                                                 }
        if( vpin->at(z).first ==  11 ){                                            vpout[vpin->at(z).second] = e.init_pz;                                                                                                 }
        if( vpin->at(z).first ==  12 ){                                            vpout[vpin->at(z).second] = e.delta_calc_energy();                                                                                     }
        if( vpin->at(z).first ==  13 ){                                            vpout[vpin->at(z).second] = e.delta_calc_prad();                                                                                       }
        if( vpin->at(z).first ==  14 ){                                            vpout[vpin->at(z).second] = e.delta_calc_plong();                                                                                      }
        if( vpin->at(z).first ==  15 ){ if( e.calc_ptot()                 != 0. ){ vpout[vpin->at(z).second] = e.delta_calc_ptot();                                                                 }  else { continue; } }
        if( vpin->at(z).first ==  16 ){                                            vpout[vpin->at(z).second] = e.delta_calc_azimuth();                                                                                    }
        if( vpin->at(z).first ==  17 ){ if( e.calc_ptot()                 != 0. ){ vpout[vpin->at(z).second] = e.delta_calc_elevation();                                                            }  else { continue; } }
        if( vpin->at(z).first ==  18 ){                                            vpout[vpin->at(z).second] = e.calc_energy();                                                                                           }
        if( vpin->at(z).first ==  19 ){                                            vpout[vpin->at(z).second] = e.calc_prad();                                                                                             }
        if( vpin->at(z).first ==  20 ){                                            vpout[vpin->at(z).second] = e.calc_plong();                                                                                            }
        if( vpin->at(z).first ==  21 ){                                            vpout[vpin->at(z).second] = e.calc_ptot();                                                                                             }
        if( vpin->at(z).first ==  22 ){                                            vpout[vpin->at(z).second] = e.calc_azimuth();                                                                                          }
        if( vpin->at(z).first ==  23 ){ if( e.calc_ptot()                 != 0. ){ vpout[vpin->at(z).second] = e.calc_elevation();                                                                  }  else { continue; } }
        if( vpin->at(z).first ==  24 ){                                            vpout[vpin->at(z).second] = e.calc_px();                                                                                               }
        if( vpin->at(z).first ==  25 ){                                            vpout[vpin->at(z).second] = e.calc_py();                                                                                               }
        if( vpin->at(z).first ==  26 ){                                            vpout[vpin->at(z).second] = e.calc_pz();                                                                                               }
        if( vpin->at(z).first ==  27 ){                                            vpout[vpin->at(z).second] = e.theo_prad();                                                                                             }
        if( vpin->at(z).first ==  28 ){                                            vpout[vpin->at(z).second] = e.theo_plong();                                                                                            }
        if( vpin->at(z).first ==  29 ){                                            vpout[vpin->at(z).second] = e.theo_azimuth();                                                                                          }
        if( vpin->at(z).first ==  30 ){ if( e.theo_ptot()                 != 0. ){ vpout[vpin->at(z).second] = e.theo_elevation();                                                                  }  else { continue; } }
        if( vpin->at(z).first ==  31 ){                                            vpout[vpin->at(z).second] = e.delta_theo_prad();                                                                                       }
        if( vpin->at(z).first ==  32 ){                                            vpout[vpin->at(z).second] = e.delta_theo_plong();                                                                                      }
        if( vpin->at(z).first ==  33 ){ if( e.theo_ptot()                 != 0. ){ vpout[vpin->at(z).second] = e.delta_theo_ptot();                                                                 }  else { continue; } }
        if( vpin->at(z).first ==  34 ){                                            vpout[vpin->at(z).second] = e.delta_theo_azimuth();                                                                                    }
        if( vpin->at(z).first ==  35 ){ if( e.theo_ptot()                 != 0. ){ vpout[vpin->at(z).second] = e.delta_theo_elevation();                                                            }  else { continue; } }
        if( vpin->at(z).first ==  36 ){ if( e.init_prad                   != 0. ){ vpout[vpin->at(z).second] = ( e.delta_calc_prad() / e.init_prad );                                               }  else { continue; } }
        if( vpin->at(z).first ==  37 ){ if( e.init_plong                  != 0. ){ vpout[vpin->at(z).second] = ( e.delta_calc_plong() / e.init_plong );                                             }  else { continue; } }
        if( vpin->at(z).first ==  38 ){ if( e.init_ptot                   != 0. ){ vpout[vpin->at(z).second] = ( e.delta_calc_ptot() / e.init_ptot );                                               }  else { continue; } }
        if( vpin->at(z).first ==  39 ){ if( e.calc_prad()                 != 0. ){ vpout[vpin->at(z).second] = ( e.delta_calc_prad() / e.calc_prad() );                                             }  else { continue; } }
        if( vpin->at(z).first ==  40 ){ if( e.calc_plong()                != 0. ){ vpout[vpin->at(z).second] = ( e.delta_calc_plong() / e.calc_plong() );                                           }  else { continue; } }
        if( vpin->at(z).first ==  41 ){ if( e.calc_ptot()                 != 0. ){ vpout[vpin->at(z).second] = ( e.delta_calc_ptot() / e.calc_ptot() );                                             }  else { continue; } }
        if( vpin->at(z).first ==  42 ){                                            vpout[vpin->at(z).second] = e.calc_prad() - e.init_prad;                                                                               }
        if( vpin->at(z).first ==  43 ){ if( e.init_prad                   != 0. ){ vpout[vpin->at(z).second] = ( e.calc_prad() - e.init_prad ) / e.init_prad;                                       }  else { continue; } }
        if( vpin->at(z).first ==  44 ){                                            vpout[vpin->at(z).second] = e.calc_plong() - e.init_plong;                                                                             }
        if( vpin->at(z).first ==  45 ){ if( e.init_plong                  != 0. ){ vpout[vpin->at(z).second] = ( e.calc_plong() - e.init_plong ) / e.init_plong;                                    }  else { continue; } }
        if( vpin->at(z).first ==  46 ){                                            vpout[vpin->at(z).second] = e.calc_ptot() - e.init_ptot;                                                                               }
        if( vpin->at(z).first ==  47 ){ if( e.init_ptot                   != 0. ){ vpout[vpin->at(z).second] = ( e.calc_ptot() - e.init_ptot ) / e.init_ptot;                                       }  else { continue; } }
        if( vpin->at(z).first ==  48 ){                                            vpout[vpin->at(z).second] = e.calc_azimuth() - e.init_azimuth();                                                                       }
        if( vpin->at(z).first ==  49 ){ if( e.calc_ptot()                 != 0. ){ vpout[vpin->at(z).second] = e.calc_elevation() - e.init_elevation();                                             }  else { continue; } }
        if( vpin->at(z).first ==  50 ){                                            vpout[vpin->at(z).second] = e.calc_energy() - e.init_ke;                                                                               }
        if( vpin->at(z).first ==  51 ){                                            vpout[vpin->at(z).second] = e.theo_prad() - e.init_prad;                                                                               }
        if( vpin->at(z).first ==  52 ){ if( e.init_prad                   != 0. ){ vpout[vpin->at(z).second] = ( e.theo_prad() - e.init_prad ) / e.init_prad;                                       }  else { continue; } }
        if( vpin->at(z).first ==  53 ){                                            vpout[vpin->at(z).second] = e.theo_plong() - e.init_plong;                                                                             }
        if( vpin->at(z).first ==  54 ){ if( e.init_plong                  != 0. ){ vpout[vpin->at(z).second] = ( e.theo_plong() - e.init_plong ) / e.init_plong;                                    }  else { continue; } }
        if( vpin->at(z).first ==  55 ){                                            vpout[vpin->at(z).second] = e.theo_ptot() - e.init_ptot;                                                                               }
        if( vpin->at(z).first ==  56 ){ if( e.init_ptot                   != 0. ){ vpout[vpin->at(z).second] = ( e.theo_ptot() - e.init_ptot ) / e.init_ptot;                                       }  else { continue; } }
        if( vpin->at(z).first ==  57 ){                                            vpout[vpin->at(z).second] = e.theo_azimuth() - e.init_azimuth();                                                                       }
        if( vpin->at(z).first ==  58 ){ if( e.theo_ptot()                 != 0. ){ vpout[vpin->at(z).second] = e.theo_elevation() - e.init_elevation();                                             }  else { continue; } }
        if( vpin->at(z).first ==  59 ){                                            vpout[vpin->at(z).second] = e.theo_energy() - e.init_ke;                                                                               }
        if( vpin->at(z).first ==  60 ){                                            vpout[vpin->at(z).second] = ( omega*e.tof - TMath::Floor( omega*e.tof/( 2*TMath::Pi() ) ) * 2*TMath::Pi() )/2.;                        }
        if( vpin->at(z).first ==  61 ){                                            vpout[vpin->at(z).second] = TMath::Sin( omega*e.tof/2. );                                                                              }
        if( vpin->at(z).first ==  62 ){                                            vpout[vpin->at(z).second] = TMath::Abs( TMath::Sin( omega*e.tof/2. ) );                                                                }
        if( vpin->at(z).first ==  63 ){                                            vpout[vpin->at(z).second] = e.splat_x;                                                                                                 }
        if( vpin->at(z).first ==  64 ){                                            vpout[vpin->at(z).second] = e.splat_y;                                                                                                 }
        if( vpin->at(z).first ==  65 ){                                            vpout[vpin->at(z).second] = e.splat_z;                                                                                                 }
        if( vpin->at(z).first ==  66 ){                                            vpout[vpin->at(z).second] = e.init_x;                                                                                                  }
        if( vpin->at(z).first ==  67 ){                                            vpout[vpin->at(z).second] = e.init_y;                                                                                                  }
        if( vpin->at(z).first ==  68 ){                                            vpout[vpin->at(z).second] = e.init_z;                                                                                                  }
        if( vpin->at(z).first ==  69 ){                                            vpout[vpin->at(z).second] = e.calc_prad() - e.theo_prad();                                                                             }
        if( vpin->at(z).first ==  70 ){                                            vpout[vpin->at(z).second] = e.calc_plong() - e.theo_plong();                                                                           }
        if( vpin->at(z).first ==  71 ){                                            vpout[vpin->at(z).second] = e.calc_ptot() - e.theo_ptot();                                                                             }
        if( vpin->at(z).first ==  72 ){                                            vpout[vpin->at(z).second] = e.calc_azimuth() - e.theo_azimuth();                                                                       }
        if( vpin->at(z).first ==  73 ){ if( e.theo_ptot() * e.calc_ptot() != 0. ){ vpout[vpin->at(z).second] = e.calc_elevation() - e.theo_elevation();                                             }  else { continue; } } // MCP MEASURES FROM HERE
        if( vpin->at(z).first ==  74 ){                                            vpout[vpin->at(z).second] = e.calc_energy() - e.theo_energy();                                                                         }
        if( vpin->at(z).first ==  75 ){                                            vpout[vpin->at(z).second] = e.delta_mcp_calc_energy();                                                                                 }
        if( vpin->at(z).first ==  76 ){                                            vpout[vpin->at(z).second] = e.delta_mcp_calc_prad();                                                                                   }
        if( vpin->at(z).first ==  77 ){                                            vpout[vpin->at(z).second] = e.delta_mcp_calc_plong();                                                                                  }
        if( vpin->at(z).first ==  78 ){ if( e.mcp_calc_ptot()             != 0. ){ vpout[vpin->at(z).second] = e.delta_mcp_calc_ptot();                                                             }  else { continue; } }
        if( vpin->at(z).first ==  79 ){                                            vpout[vpin->at(z).second] = e.delta_mcp_calc_azimuth();                                                                                }
        if( vpin->at(z).first ==  80 ){ if( e.mcp_calc_ptot()             != 0. ){ vpout[vpin->at(z).second] = e.delta_mcp_calc_elevation();                                                        }  else { continue; } }
        if( vpin->at(z).first ==  81 ){                                            vpout[vpin->at(z).second] = e.mcp_calc_energy();                                                                                       }
        if( vpin->at(z).first ==  82 ){                                            vpout[vpin->at(z).second] = e.mcp_calc_prad();                                                                                         }
        if( vpin->at(z).first ==  83 ){                                            vpout[vpin->at(z).second] = e.mcp_calc_plong();                                                                                        }
        if( vpin->at(z).first ==  84 ){                                            vpout[vpin->at(z).second] = e.mcp_calc_ptot();                                                                                         }
        if( vpin->at(z).first ==  85 ){                                            vpout[vpin->at(z).second] = e.mcp_calc_azimuth();                                                                                      }
        if( vpin->at(z).first ==  86 ){ if( e.mcp_calc_ptot()             != 0. ){ vpout[vpin->at(z).second] = e.mcp_calc_elevation();                                                              }  else { continue; } }
        if( vpin->at(z).first ==  87 ){                                            vpout[vpin->at(z).second] = e.mcp_calc_px();                                                                                           }
        if( vpin->at(z).first ==  88 ){                                            vpout[vpin->at(z).second] = e.mcp_calc_py();                                                                                           }
        if( vpin->at(z).first ==  89 ){                                            vpout[vpin->at(z).second] = e.mcp_calc_pz();                                                                                           }
        if( vpin->at(z).first ==  90 ){ if( e.init_prad                   != 0. ){ vpout[vpin->at(z).second] = ( e.delta_mcp_calc_prad() / e.init_prad );                                           }  else { continue; } }
        if( vpin->at(z).first ==  91 ){ if( e.init_plong                  != 0. ){ vpout[vpin->at(z).second] = ( e.delta_mcp_calc_plong() / e.init_plong );                                         }  else { continue; } }
        if( vpin->at(z).first ==  92 ){ if( e.init_ptot                   != 0. ){ vpout[vpin->at(z).second] = ( e.delta_mcp_calc_ptot() / e.init_ptot );                                           }  else { continue; } }
        if( vpin->at(z).first ==  93 ){ if( e.mcp_calc_prad()             != 0. ){ vpout[vpin->at(z).second] = ( e.delta_mcp_calc_prad() / e.mcp_calc_prad() );                                     }  else { continue; } }
        if( vpin->at(z).first ==  94 ){ if( e.mcp_calc_plong()            != 0. ){ vpout[vpin->at(z).second] = ( e.delta_mcp_calc_plong() / e.mcp_calc_plong() );                                   }  else { continue; } }
        if( vpin->at(z).first ==  95 ){ if( e.mcp_calc_ptot()             != 0. ){ vpout[vpin->at(z).second] = ( e.delta_mcp_calc_ptot() / e.mcp_calc_ptot() );                                     }  else { continue; } }
        if( vpin->at(z).first ==  96 ){                                            vpout[vpin->at(z).second] = e.mcp_calc_prad() - e.init_prad;                                                                           }
        if( vpin->at(z).first ==  97 ){ if( e.init_prad                   != 0. ){ vpout[vpin->at(z).second] = ( e.mcp_calc_prad() - e.init_prad ) / e.init_prad;                                   }  else { continue; } }
        if( vpin->at(z).first ==  98 ){                                            vpout[vpin->at(z).second] = e.mcp_calc_plong() - e.init_plong;                                                                         }
        if( vpin->at(z).first ==  99 ){ if( e.init_plong                  != 0. ){ vpout[vpin->at(z).second] = ( e.mcp_calc_plong() - e.init_plong ) / e.init_plong;                                }  else { continue; } }
        if( vpin->at(z).first == 100 ){                                            vpout[vpin->at(z).second] = e.mcp_calc_ptot() - e.init_ptot;                                                                           }
        if( vpin->at(z).first == 101 ){ if( e.init_ptot                   != 0. ){ vpout[vpin->at(z).second] = ( e.mcp_calc_ptot() - e.init_ptot ) / e.init_ptot;                                   }  else { continue; } }
        if( vpin->at(z).first == 102 ){                                            vpout[vpin->at(z).second] = e.mcp_calc_azimuth() - e.init_azimuth();                                                                   }
        if( vpin->at(z).first == 103 ){ if( e.mcp_calc_ptot()             != 0. ){ vpout[vpin->at(z).second] = e.mcp_calc_elevation() - e.init_elevation();                                         }  else { continue; } }
        if( vpin->at(z).first == 104 ){                                            vpout[vpin->at(z).second] = e.mcp_calc_energy() - e.init_ke;                                                                           }
        
        if( vpin->at(z).first == 105 ){                                            vpout[vpin->at(z).second] = e.vect_prad() - e.init_prad;                                                                               }
      }
      
      outx = vpout[0];
      outy = vpout[1];
      outz = vpout[2];
      outn = so->GetReadEntry();
      
      tout->Fill();
    }
  }
  
  TCanvas *outc = new TCanvas( "outc", TString::Format( "%d_%s", run_list[0], title.Data() ), 800, 800 );
  TPad *pout = new TPad( "pout", "", 0, 0, 1, 1 );
  pout->Draw();
  pout->SetRightMargin( 0.15 );
  pout->cd();
  tout->Draw( "outy:outx:outz>>hout", "", "COLZ" );
  
  hout->GetXaxis()->SetTitle( axis_x.Data() );
  hout->GetYaxis()->SetTitle( axis_y.Data() );
  hout->GetZaxis()->SetTitle( axis_z.Data() );
  hout->GetZaxis()->SetTitleOffset( 1.3 );
  hout->GetXaxis()->SetLabelSize( 0.02 );
  hout->GetYaxis()->SetLabelSize( 0.02 );
  hout->GetZaxis()->SetLabelSize( 0.02 );
  hout->SetTitle( title.Data() );
  
  if( PRODUCE_FILE ){ 
    TFile *fout = new TFile( TString::Format( "%s.root", outfile_name.Data() ), "recreate" );
    TTree *tclone = ( TTree* )tout->Clone();
    tclone->GetBranch( "outx" )->SetTitle( TString::Format( "%s/D", outputs[select_x][3].Data() ) ); tclone->GetBranch( "outx" )->SetName( outputs[select_x][3].Data() );
    tclone->GetBranch( "outy" )->SetTitle( TString::Format( "%s/D", outputs[select_y][3].Data() ) ); tclone->GetBranch( "outy" )->SetName( outputs[select_y][3].Data() );
    tclone->GetBranch( "outz" )->SetTitle( TString::Format( "%s/D", outputs[select_z][3].Data() ) ); tclone->GetBranch( "outz" )->SetName( outputs[select_z][3].Data() );
    tclone->GetLeaf  ( "outx" )->SetTitle( outputs[select_x][3].Data() )                           ; tclone->GetLeaf  ( "outx" )->SetName( outputs[select_x][3].Data() );
    tclone->GetLeaf  ( "outy" )->SetTitle( outputs[select_y][3].Data() )                           ; tclone->GetLeaf  ( "outy" )->SetName( outputs[select_y][3].Data() );
    tclone->GetLeaf  ( "outz" )->SetTitle( outputs[select_z][3].Data() )                           ; tclone->GetLeaf  ( "outz" )->SetName( outputs[select_z][3].Data() );
    tclone->SetName( TString::Format( "tout_%s", outfile_name.Data() ) );
    tclone->SetTitle( title.Data() );
    tclone->Write();
    
    TNamed *cuts_configuration = new TNamed( "cuts_configuration", cuts.Data()                 ); cuts_configuration->Write();
    TNamed *varx               = new TNamed( "varx"              , outputs[select_x][3].Data() ); varx              ->Write();
    TNamed *vary               = new TNamed( "vary"              , outputs[select_y][3].Data() ); vary              ->Write();
    TNamed *varz               = new TNamed( "varz"              , outputs[select_z][3].Data() ); varz              ->Write();
    TNamed *unitx              = new TNamed( "unitx"             , unit_x.Data()               ); unitx             ->Write();
    TNamed *unity              = new TNamed( "unity"             , unit_y.Data()               ); unity             ->Write();
    TNamed *unitz              = new TNamed( "unitz"             , unit_z.Data()               ); unitz             ->Write();
    TNamed *cutscode           = new TNamed( "cutscode"          , cutselector.c_str()         ); cutscode          ->Write();
    TNamed *axis_x_title       = new TNamed( "axis_x_title"      , axis_x.Data()               ); axis_x_title      ->Write();
    TNamed *axis_y_title       = new TNamed( "axis_y_title"      , axis_y.Data()               ); axis_y_title      ->Write();
    TNamed *axis_z_title       = new TNamed( "axis_z_title"      , axis_z.Data()               ); axis_z_title      ->Write();
    TNamed *thisto_title       = new TNamed( "thisto_title"      , title .Data()               ); thisto_title      ->Write();
    
    fout->Close();
  }
}
