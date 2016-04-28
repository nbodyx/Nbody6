
  fort.x:  List of output units for NBODY6-GPU
  --------------------------------------------

  Unit~~KZ~~~Purpose
=========================================================================
   1     1   Common save (=1: touch STOP or TCRIT; =2: main output)
   2     2   Common save (=1: every output; =2: error check & restart)
   3     3   Data bank on OUT3 (<= 2: up to NTOT, unformatted)
   4     4   Binary data analysis (suppressed); number of NS & BH
   5     5   Initial conditions (=1: Plummer model)
   6     *   Main output stream (# 6 > 1: also binaries in bodies.f)
   7     *   Zare exchange criterion (also synch.f with # 34)    
   8     8   Binary diagnostics (NEW & END)
   9     8   Full binary table on OUT9 (>= 2 without primordials)
  10    18   Hierarchical systems on HIARCH (=1/3; also # 22)
  11    23   Output of escapers on ESC (t, m, En, V, K*, NAME)
  12    19   Stellar coalescence on COAL (coal.f)
  13    19   Stellar collision on COLL (mix.f)
  14     7   Mass percentiles (# 7 > 1)
  15    19   Collision statistics (KSTAR > 12; mix.f)
  16     *   KS rectification (ksrect.f)
  17     8   Binary evolution on BINEV (# 8 > 3)
  18     8   New binaries (binout.f)
  19     9   Table of wide binaries (NB! option #9 < 0; also #8 >= 2)
  20    19   Symbiotic stars diagnostics (mdot.f)
  21    27   Diagnostics for circularizing binary (ksint.f)
  22    34   Degenerate Roche event (roche.f)
  23     8   Eigenevolution for primordial binaries (binpop.f)
  24    11   Diagnostic output of BH-stellar disruptions (#11 < 0)
  25    27   TEV look-up warning for unperturbed binary (unpert.f)
  26    19   Second coalescence channel on COAL2
  27     7   Central velocity plot (# 7 = 5; lagr.f)
  28     7   Mean mass & central distance (#7 = 6, lagr.f)  
  29    37   High-velocity special treatment (hivel.f)
  30    18   New or terminated hierarchical systems (hiarch.f)
  31     7   Mass percentiles for heavy members (# >= 2; lagr2.f)
  32     7   Mass percentiles for light members (# >= 2; lagr2.f)
  33     3   Tidal tail members (# 3 <= 3, unformatted)
  34     3   All stars + tidal tail (# 3 > 3, formatted astro units)
  35     5   Planetesimal disk analysis (# 5 = 3; adjust.f)
  36     7   Average mass in spherical shells (# 7 = 5; lagr.f)
  37     -   Not used
  38    19   Stellar evolution diagnostics (comenv.f, mix.f, trflow.f)
  39    33   Diagnostics of hardest binary (# 33 >= 2; adjust.f)
  40     *   Diagnostics of massive body in reduce.f (ARC).
  41     *   GPU overflow diagnostics (intgrt.omp.f)
  42    42   Kozai diagnostics.
  43    19   Warning of fast stellar radius increase (mdot.f)
  44    27   Rectification of induced eccentricity (# 27 = 2; impact.f)
  45    45   Plotting BH orbit(s) (NAME = 1 or 2) by BHPLOT (n = KZ(24))
  45    45   KS binary elements; also close perturber (e > 0.90, #45 > 1)
  46    45   Binary elements for close perturber (e > 0.90, #45 > 2) 
  47    13   Interstellar cloud diagnostics (output.f)
  48    18   Unstable triple check (check3.f from GPU2 adjust.f)
  49    45   Diagnostics for second innermost BH orbit (#45 > 3)
  50     -   Not used
  51     8   Reduced eccentricity to avoid collision (BINPOP; #27 = 2)
  52    14   Circular velocity of galactic halo model (#14 = 3)
  53    11   Used by NBODY7 for SPIN (ARCHAIN; #11 < 0)
  54     *   Close encounter inside 5 maximum radii (ksint.f)
  55     *   Encounter inside tidal capture radius (# 27; ksint.f)
  56     -   Plotting file for TIME, TPHYS, <R>, MASS, NCOLL (# 44 > 0)
  57     -   BH binary plotting data (ARC)
  58     -   Diagnostic ARC CHAIN info: T, IPN, NPERT, RSUM, RGRAV, GPERT
  59     -   Not used
  60     -   Not used
  61     -   Not used
  62     -   Not used
  63     -   Not used
  64     -   Not used
  65     -   Not used
  66    27   Spin overshooting (spiral.f)
  67     -   Not used
  68     -   Not used
  69     -   Not used
  70     -   Not used
  71    27   Chaos binary update skipped in hierarchy (spiral.f)
  72     -   Not used
  73     *   Diagnostics of stable hierarchies (impact.f)
  74     -   Not used
  75     *   Diagnostics for large outer eccentricity (decide.f)
  76     -   Not used
  77     *   Termination of hierarchical system (TMDIS; ksint.f)
  78    16   Modification of KS parameters (# 16 > 1; chmod.f)
  79     -   Perturbed/degenerate chain stability (cstab2.f, cstab3.f)
  80     *   Diagnostics of stable hierarchy (suppressed; impact.f)
  81     *   Triple chain stability diagnostics (chstab.f)
  82     *   Four-body stability diagnostics (cstab4.f)
  82    12   HR diagnostics of binaries (prepared in hrplot.f)
  83    12   HR diagnostics of single stars (hrplot.f)
  84     8   Binary diagnostics (# 8 > 1; bindat.f)
  85    34   Mass transfer diagnostics on ROCHE (roche.f)
  86    19   Common envelope in chain regularization (cmbody.f)
  87     8   Hierarchical data bank on HIDAT (# 8 > 3; hidat.f)
  88    48   Three-body stability comparison diagnostics
  89    37   Diagnostic output for quadruples (two types; impact.f)
  90    48   Failed three-body stability tests (> 1)
  91    19   New blue straggler formation (mix.f)
  92    34   Spin synchronization diagnostics (suppressed)
  93     4   Black hole histogram (factor of 2 from 4 Msun; #4 > 1)
  94     -   Not used
  95    27   Diagnostic spin information for WDs (spiral.f)
  96    27   Diagnostics for eccentricity decrease (spiral.f)
  97     -   Not used
  98     -   Not used
  99    22   Initial conditions in physical units (# 22 = -1/5; scale.f)
=========================================================================

SJA 28/3/14
