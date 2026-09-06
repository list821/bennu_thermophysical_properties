KPL/FK

OSIRIS-REx Target Body DSK Surface ID Codes
========================================================================

   This kernel contains a set of OSIRIS-REx Target Body DSK Surface ID
   Codes for Bennu. 


Version and Date
========================================================================

   Version 0.1 - May 7, 2019 - Michael Nolan / U Arizona

      First version. Model is "SPC v020", based on [3] and [4] and is
      delivered at ~3-m sample spacing.

   Version 0.2 - Aug 5, 2020 - Michael Nolan / U Arizona and Ray 
      Espiritu / Johns Hopkin University Applied Physics Laboratory.

      Added model "ALT v020", based on [5], and is delivered at 3-m and
      80-cm spacing globally and 5-cm near the planned TAG sites.

   Version 0.3 - Jul 9, 2021 - Michael Nolan / U Arizona and Ray 
      Espiritu / Johns Hopkin University Applied Physics Laboratory.

      Added model "ALT v021", based on [5], and is delivered at 12.5-m,
      6m, 3m, 1.5m, and 80cm spacing globally and 5cm near the
      planned TAG sites.
      
      In addition there is a 40cm Poisson reconstructed global model and
      5cm Poisson reconstructed models near the planned TAG sites.

References
========================================================================

   1. ``Kernel Pool Required Reading''

   2. ``DS-Kernel Required Reading''

   3. SPC Shape Model Evaluation 2019/01/21 (C. Johnson, L. Philpott,
      M. Al Asad, O. Barnouin)
 
   4. O.S. Barnouni et al. (2019) Shape of (101955) Bennu indicitive of
      a rubble pile with internal stiffness. Nature Geoscience 12,
      247-252

   5. M.G. Daly et al. (2020) Hemispherical Differences in the Shape and
      Topography of Asteroid (101955) Bennu. Science Advances, 6,
      eabd3649.


Contact Information
========================================================================

      NAIF at JPL:

           Boris Semenov
           (818) 354-8136
           Boris.Semenov@jpl.nasa.gov


Implementation Notes
========================================================================

   This file is used by the SPICE system as follows: programs that make
   use of this frame kernel must `load' the kernel, normally during
   program initialization. The SPICELIB routine FURNSH and CSPICE
   function furnsh_c load a kernel file into the kernel pool as shown
   below.

      CALL FURNSH ( 'frame_kernel_name' )
      furnsh_c    ( "frame_kernel_name" );

   This file was created and may be updated with a text editor or word
   processor.


Definition Section
========================================================================

   This section contains name to ID mappings for the OSIRIS-REx target
   body DSK surfaces. These mappings are supported by all SPICE
   toolkits with integrated DSK capabilities (version N0066 or later).

   The IDs have the form PRRRRRVVV, where P is a 1 digit kind / producer
   / technique identifier (tru=1, spc=2, ola=3), RRRRR is a 5 digit
   resolution in mm (1254, 0790, 0030, 0005), VVV is a 3 digit version
   (120, 450). Note that each producer may have models with the same
   version numbers.

Version 0.1

   The previous Version 0.1 contains the project-official shape
   model "SPC v020" generated Jan 21 2019. That model was produced at
   multiple resolutions; this version is the model at ~ 3m ground
   sample distance. The model is described in [3] and [4].

Version 0.2

   Version 0.2 contains the contents of version 0.1 and
   adds project-official shape model "ALT v020" generated April 14 2020.
   That model was produced at multiple resolutions using OLA laser
   altimetry as described in [5]. These versions are global models at
   ~3m and ~88-cm spacing and local models near the selected TAG sites
   at 5-cm spacing. Note that the SPCv020 model and the ALTv020 model
   are completely different: The version numbers being the same is a
   coincidence.

Version 0.3

   This current version 0.3 contains the contents of version 0.2 and
   adds project-official shape model "ALT v021" generated Jul 9 2021.
   That model was produced at multiple resolutions using OLA laser
   altimetry as described in [5]. These versions are global models at
   ~12.5m, ~6m, ~3m, ~1.5m, and ~88cm spacing and local models near 
   the selected TAG sites at 5cm spacing.     
      
   In addition there are models that were created using a Poisson 
   reconstruction method: a 40cm global model and local models near
   the selected TAG sites at 5cm spacing.


   Asteroid Bennu Surface name/IDs:

          DSK Surface Name           ID       Body ID
      ==========================   =========  =======

      ORX_BENNU_03170MM_SPC_V020   203170020  2101955
      ORX_BENNU_03250MM_ALT_V020   303250020  2101955
      ORX_BENNU_00880MM_ALT_V020   300880020  2101955
      ORX_BENNU_00050MM_ALT_V020   300050020  2101955

   Name-ID Mapping keywords:

   \begindata

      NAIF_SURFACE_NAME += 'ORX_BENNU_03170MM_SPC_V020'
      NAIF_SURFACE_CODE += 203170020
      NAIF_SURFACE_BODY += 2101955

      NAIF_SURFACE_NAME += 'ORX_BENNU_03250MM_ALT_V020'
      NAIF_SURFACE_CODE += 303250020
      NAIF_SURFACE_BODY += 2101955

      NAIF_SURFACE_NAME += 'ORX_BENNU_00880MM_ALT_V020'
      NAIF_SURFACE_CODE += 300880020
      NAIF_SURFACE_BODY += 2101955

      NAIF_SURFACE_NAME += 'ORX_BENNU_00050MM_ALT_V020'
      NAIF_SURFACE_CODE += 300050020
      NAIF_SURFACE_BODY += 2101955

   \begintext

   
   Asteroid Bennu Surface name/IDs:

          DSK Surface Name           ID       Body ID
      ==========================   =========  =======

      ORX_BENNU_12600MM_ALT_V021   312600021  2101955
      ORX_BENNU_06370MM_ALT_V021   306370021  2101955
      ORX_BENNU_03250MM_ALT_V021   303250021  2101955
      ORX_BENNU_01680MM_ALT_V021   301680021  2101955
      ORX_BENNU_00880MM_ALT_V021   300880021  2101955
      ORX_BENNU_00400MM_ALT_V021   300400021  2101955
      ORX_BENNU_00050MM_ALT_V021   300050021  2101955

   Name-ID Mapping keywords:

   \begindata

      NAIF_SURFACE_NAME += 'ORX_BENNU_12600MM_ALT_V021'
      NAIF_SURFACE_CODE += 312600021
      NAIF_SURFACE_BODY += 2101955

      NAIF_SURFACE_NAME += 'ORX_BENNU_06370MM_ALT_V021'
      NAIF_SURFACE_CODE += 306370021
      NAIF_SURFACE_BODY += 2101955

      NAIF_SURFACE_NAME += 'ORX_BENNU_03250MM_ALT_V021'
      NAIF_SURFACE_CODE += 303250021
      NAIF_SURFACE_BODY += 2101955

      NAIF_SURFACE_NAME += 'ORX_BENNU_01680MM_ALT_V021'
      NAIF_SURFACE_CODE += 301680021
      NAIF_SURFACE_BODY += 2101955

      NAIF_SURFACE_NAME += 'ORX_BENNU_00880MM_ALT_V021'
      NAIF_SURFACE_CODE += 300880021
      NAIF_SURFACE_BODY += 2101955

      NAIF_SURFACE_NAME += 'ORX_BENNU_00400MM_ALT_V021'
      NAIF_SURFACE_CODE += 300400021
      NAIF_SURFACE_BODY += 2101955

      NAIF_SURFACE_NAME += 'ORX_BENNU_00050MM_ALT_V021'
      NAIF_SURFACE_CODE += 300050021
      NAIF_SURFACE_BODY += 2101955

   \begintext


End of FK.
