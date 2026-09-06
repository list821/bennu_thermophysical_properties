KPL/MK

   This meta-kernel lists the OSIRIS-REx SPICE kernels, including OLA
   CKs, providing coverage for 2018. All of the kernels listed below
   are archived in the PDS OSIRIS-REx SPICE kernel archive. This set of
   files and the order in which they are listed were picked to provide
   the best available data and the most complete coverage for the
   specified year based on the information about the kernels available
   at the time this meta-kernel was made. For detailed information
   about the kernels listed below refer to the internal comments
   included in the kernels and the documentation accompanying the
   OSIRIS-REx SPICE kernel archive.

   It is recommended that users make a local copy of this file and
   modify the value of the PATH_VALUES keyword to point to the actual
   location of the OSIRIS-REx SPICE kernel archives' ``spice_kernels''
   directory on their system. Replacing ``/'' with ``\'' and converting
   line terminators to the format native to the user's system may also
   be required if this meta-kernel is to be used on a non-UNIX
   workstation.

   This file was created on November 6, 2023 by Alyssa Bailey, NAIF/JPL.
   The original name of this file was orx_2018_v09.tm.

   \begindata

      PATH_VALUES     = ( '..'      )

      PATH_SYMBOLS    = ( 'KERNELS' )

      KERNELS_TO_LOAD = (

                          '$KERNELS/lsk/naif0012.tls'

                          '$KERNELS/pck/pck00010.tpc'
                          '$KERNELS/pck/bennu_v17.tpc'

                          '$KERNELS/fk/orx_v14.tf'
                          '$KERNELS/fk/orx_shape_v03.tf'

                          '$KERNELS/ik/orx_lidar_v00.ti'
                          '$KERNELS/ik/orx_navcam_v02.ti'
                          '$KERNELS/ik/orx_ocams_v07.ti'
                          '$KERNELS/ik/orx_ola_v01.ti'
                          '$KERNELS/ik/orx_otes_v00.ti'
                          '$KERNELS/ik/orx_ovirs_v00.ti'
                          '$KERNELS/ik/orx_rexis_v01.ti'
                          '$KERNELS/ik/orx_stowcam_v00.ti'
                          '$KERNELS/ik/orx_struct_v00.ti'

                          '$KERNELS/sclk/orx_sclkscet_00093.tsc'

                          '$KERNELS/spk/de424.bsp'

                          '$KERNELS/spk/bennu_refdrmc_v1.bsp'

                          '$KERNELS/spk/orx_struct_v04.bsp'

                          '$KERNELS/spk/orx_170923_180710_180321_od031_v1.bsp'
                          '$KERNELS/spk/orx_180301_181201_180921_od044_v1.bsp'
                          '$KERNELS/spk/orx_180801_190302_181218_od077_v1.bsp'
                          '$KERNELS/spk/orx_181203_190302_190104_od085_v1.bsp'

                          '$KERNELS/ck/orx_ola_180125_scil2id00009.bc'
                          '$KERNELS/ck/orx_ola_180125_scil2id00010.bc'
                          '$KERNELS/ck/orx_ola_180125_scil2id00011.bc'
                          '$KERNELS/ck/orx_ola_180125_scil2id00013.bc'
                          '$KERNELS/ck/orx_ola_180312_scil2id00015.bc'
                          '$KERNELS/ck/orx_ola_180312_scil2id00016.bc'
                          '$KERNELS/ck/orx_ola_180312_scil2id00017.bc'
                          '$KERNELS/ck/orx_ola_180312_scil2id00019.bc'
                          '$KERNELS/ck/orx_ola_180312_scil2id00021.bc'
                          '$KERNELS/ck/orx_ola_180312_scil2id00022.bc'
                          '$KERNELS/ck/orx_ola_180312_scil2id00023.bc'
                          '$KERNELS/ck/orx_ola_180312_scil2id00025.bc'
                          '$KERNELS/ck/orx_ola_180716_scil2id00030.bc'
                          '$KERNELS/ck/orx_ola_180716_scil2id00031.bc'
                          '$KERNELS/ck/orx_ola_180716_scil2id00032.bc'
                          '$KERNELS/ck/orx_ola_180716_scil2id00034.bc'
                          '$KERNELS/ck/orx_ola_180717_scil2id00036.bc'
                          '$KERNELS/ck/orx_ola_181204_scil2id01000_v02.bc'
                          '$KERNELS/ck/orx_ola_181204_scil2id01001_v02.bc'
                          '$KERNELS/ck/orx_ola_181208_scil2id01002_v02.bc'
                          '$KERNELS/ck/orx_ola_181208_scil2id01003_v02.bc'
                          '$KERNELS/ck/orx_ola_181212_scil2id01004_v02.bc'
                          '$KERNELS/ck/orx_ola_181212_scil2id01005_v02.bc'
                          '$KERNELS/ck/orx_ola_181216_scil2id01006_v02.bc'
                          '$KERNELS/ck/orx_ola_181216_scil2id01007_v02.bc'

                          '$KERNELS/ck/orx_sa_rel_180101_180107_v02.bc'
                          '$KERNELS/ck/orx_sa_rel_180108_180114_v01.bc'
                          '$KERNELS/ck/orx_sa_rel_180115_180121_v01.bc'
                          '$KERNELS/ck/orx_sa_rel_180122_180128_v01.bc'
                          '$KERNELS/ck/orx_sa_rel_180129_180204_v01.bc'
                          '$KERNELS/ck/orx_sa_rel_180205_180211_v01.bc'
                          '$KERNELS/ck/orx_sa_rel_180212_180218_v01.bc'
                          '$KERNELS/ck/orx_sa_rel_180219_180225_v01.bc'
                          '$KERNELS/ck/orx_sa_rel_180226_180304_v01.bc'
                          '$KERNELS/ck/orx_sa_rel_180305_180311_v02.bc'
                          '$KERNELS/ck/orx_sa_rel_180312_180318_v02.bc'
                          '$KERNELS/ck/orx_sa_rel_180319_180325_v01.bc'
                          '$KERNELS/ck/orx_sa_rel_180326_180401_v01.bc'
                          '$KERNELS/ck/orx_sa_rel_180402_180408_v01.bc'
                          '$KERNELS/ck/orx_sa_rel_180409_180415_v01.bc'
                          '$KERNELS/ck/orx_sa_rel_180416_180422_v01.bc'
                          '$KERNELS/ck/orx_sa_rel_180423_180429_v02.bc'
                          '$KERNELS/ck/orx_sa_rel_180430_180506_v02.bc'
                          '$KERNELS/ck/orx_sa_rel_180507_180513_v02.bc'
                          '$KERNELS/ck/orx_sa_rel_180514_180520_v02.bc'
                          '$KERNELS/ck/orx_sa_rel_180521_180527_v01.bc'
                          '$KERNELS/ck/orx_sa_rel_180528_180603_v02.bc'
                          '$KERNELS/ck/orx_sa_rel_180604_180610_v02.bc'
                          '$KERNELS/ck/orx_sa_rel_180611_180617_v01.bc'
                          '$KERNELS/ck/orx_sa_rel_180618_180624_v02.bc'
                          '$KERNELS/ck/orx_sa_rel_180625_180701_v02.bc'
                          '$KERNELS/ck/orx_sa_rel_180702_180708_v02.bc'
                          '$KERNELS/ck/orx_sa_rel_180709_180715_v02.bc'
                          '$KERNELS/ck/orx_sa_rel_180716_180722_v02.bc'
                          '$KERNELS/ck/orx_sa_rel_180723_180729_v01.bc'
                          '$KERNELS/ck/orx_sa_rel_180730_180805_v01.bc'
                          '$KERNELS/ck/orx_sa_rel_180806_180812_v02.bc'
                          '$KERNELS/ck/orx_sa_rel_180813_180819_v02.bc'
                          '$KERNELS/ck/orx_sa_rel_180820_180826_v01.bc'
                          '$KERNELS/ck/orx_sa_rel_180827_180902_v01.bc'
                          '$KERNELS/ck/orx_sa_rel_180903_180909_v02.bc'
                          '$KERNELS/ck/orx_sa_rel_180910_180916_v01.bc'
                          '$KERNELS/ck/orx_sa_rel_180917_180923_v02.bc'
                          '$KERNELS/ck/orx_sa_rel_180924_180930_v02.bc'
                          '$KERNELS/ck/orx_sa_rel_181001_181007_v01.bc'
                          '$KERNELS/ck/orx_sa_rel_181008_181014_v02.bc'
                          '$KERNELS/ck/orx_sa_rel_181015_181021_v02.bc'
                          '$KERNELS/ck/orx_sa_rel_181022_181028_v02.bc'
                          '$KERNELS/ck/orx_sa_rel_181029_181104_v02.bc'
                          '$KERNELS/ck/orx_sa_rel_181105_181111_v02.bc'
                          '$KERNELS/ck/orx_sa_rel_181112_181118_v01.bc'
                          '$KERNELS/ck/orx_sa_rel_181119_181125_v01.bc'
                          '$KERNELS/ck/orx_sa_rel_181126_181202_v01.bc'
                          '$KERNELS/ck/orx_sa_rel_181203_181209_v01.bc'
                          '$KERNELS/ck/orx_sa_rel_181210_181216_v01.bc'
                          '$KERNELS/ck/orx_sa_rel_181217_181223_v01.bc'
                          '$KERNELS/ck/orx_sa_rel_181224_181230_v01.bc'
                          '$KERNELS/ck/orx_sa_rel_181231_190106_v01.bc'

                          '$KERNELS/ck/orx_sc_rel_180101_180107_v02.bc'
                          '$KERNELS/ck/orx_sc_rel_180108_180114_v01.bc'
                          '$KERNELS/ck/orx_sc_rel_180115_180121_v01.bc'
                          '$KERNELS/ck/orx_sc_rel_180122_180128_v01.bc'
                          '$KERNELS/ck/orx_sc_rel_180129_180204_v01.bc'
                          '$KERNELS/ck/orx_sc_rel_180205_180211_v01.bc'
                          '$KERNELS/ck/orx_sc_rel_180212_180218_v01.bc'
                          '$KERNELS/ck/orx_sc_rel_180219_180225_v01.bc'
                          '$KERNELS/ck/orx_sc_rel_180226_180304_v01.bc'
                          '$KERNELS/ck/orx_sc_rel_180305_180311_v02.bc'
                          '$KERNELS/ck/orx_sc_rel_180312_180318_v02.bc'
                          '$KERNELS/ck/orx_sc_rel_180319_180325_v01.bc'
                          '$KERNELS/ck/orx_sc_rel_180326_180401_v01.bc'
                          '$KERNELS/ck/orx_sc_rel_180402_180408_v01.bc'
                          '$KERNELS/ck/orx_sc_rel_180409_180415_v01.bc'
                          '$KERNELS/ck/orx_sc_rel_180416_180422_v01.bc'
                          '$KERNELS/ck/orx_sc_rel_180423_180429_v02.bc'
                          '$KERNELS/ck/orx_sc_rel_180430_180506_v02.bc'
                          '$KERNELS/ck/orx_sc_rel_180507_180513_v02.bc'
                          '$KERNELS/ck/orx_sc_rel_180514_180520_v02.bc'
                          '$KERNELS/ck/orx_sc_rel_180521_180527_v01.bc'
                          '$KERNELS/ck/orx_sc_rel_180528_180603_v02.bc'
                          '$KERNELS/ck/orx_sc_rel_180604_180610_v02.bc'
                          '$KERNELS/ck/orx_sc_rel_180611_180617_v01.bc'
                          '$KERNELS/ck/orx_sc_rel_180618_180624_v02.bc'
                          '$KERNELS/ck/orx_sc_rel_180625_180701_v02.bc'
                          '$KERNELS/ck/orx_sc_rel_180702_180708_v02.bc'
                          '$KERNELS/ck/orx_sc_rel_180709_180715_v02.bc'
                          '$KERNELS/ck/orx_sc_rel_180716_180722_v02.bc'
                          '$KERNELS/ck/orx_sc_rel_180723_180729_v01.bc'
                          '$KERNELS/ck/orx_sc_rel_180730_180805_v01.bc'
                          '$KERNELS/ck/orx_sc_rel_180806_180812_v02.bc'
                          '$KERNELS/ck/orx_sc_rel_180813_180819_v02.bc'
                          '$KERNELS/ck/orx_sc_rel_180820_180826_v01.bc'
                          '$KERNELS/ck/orx_sc_rel_180827_180902_v01.bc'
                          '$KERNELS/ck/orx_sc_rel_180903_180909_v02.bc'
                          '$KERNELS/ck/orx_sc_rel_180910_180916_v01.bc'
                          '$KERNELS/ck/orx_sc_rel_180917_180923_v02.bc'
                          '$KERNELS/ck/orx_sc_rel_180924_180930_v02.bc'
                          '$KERNELS/ck/orx_sc_rel_181001_181007_v02.bc'
                          '$KERNELS/ck/orx_sc_rel_181008_181014_v02.bc'
                          '$KERNELS/ck/orx_sc_rel_181015_181021_v02.bc'
                          '$KERNELS/ck/orx_sc_rel_181022_181028_v02.bc'
                          '$KERNELS/ck/orx_sc_rel_181029_181104_v02.bc'
                          '$KERNELS/ck/orx_sc_rel_181105_181111_v02.bc'
                          '$KERNELS/ck/orx_sc_rel_181112_181118_v01.bc'
                          '$KERNELS/ck/orx_sc_rel_181119_181125_v01.bc'
                          '$KERNELS/ck/orx_sc_rel_181126_181202_v01.bc'
                          '$KERNELS/ck/orx_sc_rel_181203_181209_v01.bc'
                          '$KERNELS/ck/orx_sc_rel_181210_181216_v01.bc'
                          '$KERNELS/ck/orx_sc_rel_181217_181223_v01.bc'
                          '$KERNELS/ck/orx_sc_rel_181224_181230_v01.bc'
                          '$KERNELS/ck/orx_sc_rel_181231_190106_v01.bc'

                          '$KERNELS/ck/orx_struct_mapcam_v01.bc'
                          '$KERNELS/ck/orx_struct_polycam_v01.bc'

             '$KERNELS/dsk/bennu_g_00880mm_alt_obj_0000n00000_v021a.bds'

                        )

   \begintext

End of MK file.
