J/A+A/662/A40           Hot subdwarf stars studied with Gaia     (Culpan+, 2022)
The population of hot subdwarf stars studied with Gaia.
IV. Catalogues of hot subluminous stars based on Gaia EDR3.
    Culpan R., Geier S., Reindl N., Pelisoli I., Gentile Fusillo N.,
    Vorontseva A.
    <Astron. Astrophys. 662, A40 (2022)>
    =2022A&A...662A..40C        (SIMBAD/NED BibCode)
ADC_Keywords: Stars, subdwarf ; Photometry
Keywords: subdwarfs - Hertzsprung-Russel and C-M diagrams - binaries: general -
          stars: horizontal branch - catalogs

Abstract:
    In light of substantial new discoveries of hot subdwarfs by ongoing
    spectroscopic surveys and the availability of the Gaia mission Early
    Data Release 3 (EDR3), we compiled new releases of two catalogues of
    hot subluminous stars: The data release 3 (DR3) catalogue of the known
    hot subdwarf stars contains 6616 unique sources and provides
    multi-band photometry, and astrometry from Gaia EDR3 as well as
    classifications based on spectroscopy and colours. This is an increase
    of 742 objects over the DR2 catalogue. This new catalogue provides
    atmospheric parameters of 3087 stars and radial velocities of 2791
    stars from the literature. In addition, we have updated the Gaia Data
    Release 2 (DR2) catalogue of hot subluminous stars using the improved
    accuracy of the Gaia EDR3 data set together with updated quality and
    selection criteria to produce the Gaia EDR3 catalogue of 61,585 hot
    subluminous stars, representing an increase of 21785 objects. The
    improvements, compared to Gaia DR2, in Gaia EDR3's astrometry and
    photometry have enabled us to define more sophisticated selection
    functions. In particular, we improved hot subluminous star detection
    in the crowded regions of the Galactic Plane as well as in the
    direction of the Magellanic Clouds by including sources with close
    apparent neighbours but with flux levels that dominate the
    neighbourhood.

Description:
    Here we present new releases of the catalogue of known hot subdwarfs
    (DR3) and the catalogue of hot subluminous star candidates (DR2) based
    on improved data from Gaia EDR3 (Gaia Collaboration,
    2021A&A...649A...1G, Cat. I.350), new results from spectroscopic
    surveys, and an extensive literature search.

File Summary:
--------------------------------------------------------------------------------
 FileName      Lrecl  Records   Explanations
--------------------------------------------------------------------------------
ReadMe            80        .   This file
hotsd.dat        383    61585   Hot SD Gaia EDR3 catalogue
knownhsd.dat    1012     6616   Known hotSD catalogue
--------------------------------------------------------------------------------

See also:
        I/350 : Gaia EDR3 (Gaia Collaboration, 2020)
 J/A+A/649/A3 : Gaia Early Data Release 3 photometric passbands (Riello+, 2021)

Byte-by-byte Description of file: hotsd.dat
--------------------------------------------------------------------------------
   Bytes Format Units     Label     Explanations
--------------------------------------------------------------------------------
   1- 19  I19   ---       GaiaEDR3  Gaia EDR3 source_id
  21- 42 F22.18 deg       RAdeg     Right ascension (ICRS) at Ep=2016.0
  44- 65 E22.18 deg       DEdeg     Declination (ICRS) at Ep=2016.0
  67- 89 F23.19 deg       GLON      Galactic longitude
  91-112 F22.18 deg       GLAT      Galactic latitude
 114-135 E22.18 mas       Plx       Gaia EDR3 parallax
 137-148 F12.10 mas     e_Plx       Error on Gaia parallax
 150-169 F20.17 mag       GMAG      ?=- Absolute G magnitude
 171-180  F10.7 mag       Gmag      Gaia apparent G magnitude
 182-195 E14.10 mag       BP-RP     Gaia G(BP)-G(RP) colour
 197-218 E22.18 ---       E(BP/RP)c ?=- Gaia colour excess corrected
                                     (Riello et al., 2021A&A...649A...3R,
                                     Cat. J/A+A/649/A3)
 220-241 E22.18 mas/yr    pmRA      Proper motion along RA vector, pmRA*cosDE
 243-254 F12.10 mas/yr  e_pmRA      Error in proper motion along RA vector
 256-277 E22.18 mas/yr    pmDE      Proper motion along DE vector
 279-290 F12.10 mas/yr  e_pmDE      Error in proper motion along DE vector
 292-304  F13.9 mas/yr    pm        Total proper motion
 306-324 F19.16 mas/yr  e_pm        Error in total proper motion
 326-336  F11.8 ---       RUWE      Renormalised unit weight error
 338-357 F20.17 ---       Rpm       Reduced proper motion
 359-379 E21.17 ---     e_EFlux     ?=- Excess flux error
     381  I1    ---     f_Plx       [0/1] Parallax selection flag
     383  I1    ---     f_pm        [0/1] Proper motion selection flag
--------------------------------------------------------------------------------

Byte-by-byte Description of file: knownhsd.dat
--------------------------------------------------------------------------------
   Bytes   Format Units     Label      Explanations
--------------------------------------------------------------------------------
    1-  27  A27   ---       Name       Name
   29-  37  A9    ---       ---        [Gaia EDR3 -]
   39-  57  I19   ---       GaiaEDR3   ? Gaia EDR3 source_id
   59-  73 F15.11 deg       RAdeg      Right ascension (ICRS) at Ep=2016.0
   75-  89 F15.11 deg       DEdeg      Declination (ICRS) at Ep=2016.0
   91- 105 F15.11 deg       GLON       ?=- Galactic longitude
  107- 121 F15.11 deg       GLAT       ?=- Galactic latitude
  123- 135  A13   ---       SpClass    Spectral classification
  137- 154  A18   ---       SpClassS   Simbad spectral classification
  156- 160  A5    ---       CSDSS      SDSS colours
  162- 166  A5    ---       CAPASS     APASS colours
  168- 172  A5    ---       CPS1       PS1 colours
  174- 178  A5    ---       CSKYM      SKYM colours
  180- 186  E7.4  mas       Plx        ?=- Parallax
  188- 194  E7.4  mas       PlxZP      ?=- Parallax ZP
  196- 201  F6.4  mas     e_Plx        ?=- Parallax error
  203- 209  F7.4  mag       GMAG       ?=- Absolute Gaia G magnitude
  211- 217  F7.4  mag       GGAIA      ?=- Gaia G magnitude
  219- 224  F6.4  mag     e_GGAIA      ?=- Gaia G magnitude error
  226- 232  F7.4  mag       BPGAIA     ?=- Gaia BP magnitude
  234- 239  F6.4  mag     e_BPGAIA     ?=- Gaia BP magnitude error
  241- 247  F7.4  mag       RPGAIA     ?=- Gaia RP magnitude
  249- 254  F6.4  mag     e_RPGAIA     ?=- Gaia RP magnitude error
  256- 263  F8.4  mas/yr    pmRAGAIA   ?=- Gaia proper motion along RA,pmRA*cosDE
  265- 270  F6.4  mas/yr  e_pmRAGAIA   ?=- Gaia proper motion along RA error
  272- 280  F9.4  mas/yr    pmDEGAIA   ?=- Gaia proper motion along DE
  282- 287  F6.4  mas/yr  e_pmDEGAIA   ?=- Gaia proper motion along DE error
  289- 292  I4    km/s      RVSDSS     ?=- SDSS radial velocity
  294- 295  I2    km/s    e_RVSDSS     ?=- SDSS radial velocity error
  297- 300  I4    km/s      RVLAMOST   ?=- LAMOS radial velocity
  302- 304  I3    km/s    e_RVLAMOST   ?=- LAMOS radial velocity error
  306- 311  I6    K         Teff       ?=- Effective temperature
  313- 317  I5    K       e_Teff       ?=- Effective temperature error
  319- 322  F4.2  [cm/s2]   logg       ?=- Surface gravity
  324- 327  F4.2  [cm/s2] e_logg       ?=- Surface gravity error
  329- 333  F5.2  [-]       logY       ?=- Stellar surface He abundance,
                                        Y=n(He)/n(H)
  335- 339  F5.2  [-]     e_logY       []?=- Stellar surface He abundance error
  341- 383  A43   ---       Ref        Parameters references
  385- 390  F6.3  mag       E(B-V)     Interstellar redenning
  392- 396  F5.3  mag     e_E(B-V)     Interstellar redenning error
  398- 404  F7.3  mag       AV         Absorption in V band
  406- 411  F6.3  mag       FUVGALEX   ?=- FUV GALEX magnitude
  413- 417  F5.3  mag     e_FUVGALEX   ?=- FUV GALEX magnitude error
  419- 424  F6.3  mag       NUVGALEX   ?=- NUV GALEX magnitude
  426- 430  F5.3  mag     e_NUVGALEX   ?=- NUV GALEX magnitude error
  432- 437  F6.3  mag       VAPASS     ?=- V APASS magnitude
  439- 443  F5.3  mag     e_VAPASS     ?=- V APASS magnitude error
  445- 450  F6.3  mag       BAPASS     ?=- B APASS magnitude
  452- 456  F5.3  mag     e_BAPASS     ?=- B APASS magnitude error
  458- 463  F6.3  mag       gAPASS     ?=- g APASS magnitude
  465- 469  F5.3  mag     e_gAPASS     ?=- g APASS magnitude error
  471- 476  F6.3  mag       rAPASS     ?=- r APASS magnitude
  478- 482  F5.3  mag     e_rAPASS     ?=- r APASS magnitude error
  484- 489  F6.3  mag       iAPASS     ?=- i APASS magnitude
  491- 495  F5.3  mag     e_iAPASS     ?=- i APASS magnitude error
  497- 502  F6.3  mag       uSDSS      ?=- u SDSS magnitude
  504- 509  F6.3  mag     e_uSDSS      ?=- u SDSS magnitude error
  511- 516  F6.3  mag       gSDSS      ?=- g SDSS magnitude
  518- 522  F5.3  mag     e_gSDSS      ?=- g SDSS magnitude error
  524- 529  F6.3  mag       rSDSS      ?=- r SDSS magnitude
  531- 535  F5.3  mag     e_rSDSS      ?=- r SDSS magnitude error
  537- 542  F6.3  mag       iSDSS      ?=- i SDSS magnitude
  544- 548  F5.3  mag     e_iSDSS      ?=- i SDSS magnitude error
  550- 555  F6.3  mag       zSDSS      ?=- z SDSS magnitude
  557- 561  F5.3  mag     e_zSDSS      ?=- z SDSS magnitude error
  563- 569  E7.4  mag       uVST       ?=- u VST magnitude
  571- 577  E7.2  mag     e_uVST       ?=- u VST magnitude error
  579- 585  F7.4  mag       gVST       ?=- g VST magnitude
  587- 593  E7.2  mag     e_gVST       ?=- g VST magnitude error
  595- 601  F7.4  mag       rVST       ?=- r VST magnitude
  603- 609  E7.2  mag     e_rVST       ?=- r VST magnitude error
  611- 617  F7.4  mag       iVST       ?=- i VST magnitude
  619- 625  E7.2  mag     e_iVST       ?=- i VST magnitude error
  627- 633  F7.4  mag       zVST       ?=- z VST magnitude
  635- 641  E7.5  mag     e_zVST       ?=- z VST magnitude error
  643- 648  F6.3  mag       uSKYM      ?=- u SKYM magnitude
  650- 654  F5.3  mag     e_uSKYM      ?=- u SKYM magnitude error
  656- 661  F6.3  mag       vSKYM      ?=- v SKYM magnitude
  663- 667  F5.3  mag     e_vSKYM      ?=- v SKYM magnitude error
  669- 674  F6.3  mag       gSKYM      ?=- g SKYM magnitude
  676- 680  F5.3  mag     e_gSKYM      ?=- g SKYM magnitude error
  682- 687  F6.3  mag       rSKYM      ?=- r SKYM magnitude
  689- 693  F5.3  mag     e_rSKYM      ?=- r SKYM magnitude error
  695- 700  F6.3  mag       iSKYM      ?=- i SKYM magnitude
  702- 706  F5.3  mag     e_iSKYM      ?=- i SKYM magnitude error
  708- 713  F6.3  mag       zSKYM      ?=- z SKYM magnitude
  715- 719  F5.3  mag     e_zSKYM      ?=- z SKYM magnitude error
  721- 727  F7.4  mag       gPS1       ?=- g PS1 magnitude
  729- 734  E6.4  mag     e_gPS1       ?=- g PS1 magnitude error
  736- 742  F7.4  mag       rPS1       ?=- r PS1 magnitude
  744- 750  E7.4  mag     e_rPS1       ?=- r PS1 magnitude error
  751- 757  F7.4  mag       iPS1       ?=- i PS1 magnitude
  759- 764  E6.4  mag     e_iPS1       ?=- i PS1 magnitude error
  766- 772  F7.4  mag       zPS1       ?=- z PS1 magnitude
  774- 780  E7.4  mag     e_zPS1       ?=- z PS1 magnitude error
  781- 787  F7.4  mag       yPS1       ?=- y PS1 magnitude
  789- 794  E6.4  mag     e_yPS1       ?=- y PS1 magnitude error
  796- 801  F6.3  mag       J2MASS     ?=- J 2MASS magnitude
  803- 807  F5.3  mag     e_J2MASS     ?=- J 2MASS magnitude error
  809- 814  F6.3  mag       H2MASS     ?=- H 2MASS magnitude
  816- 820  F5.3  mag     e_H2MASS     ?=- H 2MASS magnitude error
  822- 827  F6.3  mag       K2MASS     ?=- K 2MASS magnitude
  829- 833  F5.3  mag     e_K2MASS     ?=- K 2MASS magnitude error
  835- 840  F6.3  mag       YUKIDSS    ?=- Y UKIDSS magnitude
  842- 846  F5.3  mag     e_YUKIDSS    ?=- Y UKIDSS magnitude error
  848- 853  F6.3  mag       JUKIDSS    ?=- J UKIDSS magnitude
  855- 859  F5.3  mag     e_JUKIDSS    ?=- J UKIDSS magnitude error
  861- 866  F6.3  mag       HUKIDSS    ?=- H UKIDSS magnitude
  868- 872  F5.3  mag     e_HUKIDSS    ?=- H UKIDSS magnitude error
  874- 879  F6.3  mag       KUKIDSS    ?=- K UKIDSS magnitude
  881- 885  F5.3  mag     e_KUKIDSS    ?=- K UKIDSS magnitude error
  887- 893  F7.4  mag       ZVISTA     ?=- Z VISTA magnitude
  895- 900  E6.4  mag     e_ZVISTA     ?=- Z VISTA magnitude error
  902- 908  F7.4  mag       YVISTA     ?=- Y VISTA magnitude
  910- 915  E6.4  mag     e_YVISTA     ?=- Y VISTA magnitude error
  917- 923  F7.4  mag       JVISTA     ?=- J VISTA magnitude
  925- 930  E6.4  mag     e_JVISTA     ?=- J VISTA magnitude error
  932- 938  F7.4  mag       HVISTA     ?=- H VISTA magnitude
  940- 945  E6.4  mag     e_HVISTA     ?=- H VISTA magnitude error
  947- 953  F7.4  mag       KsVISTA    ?=- Ks VISTA magnitude
  955- 960  E6.2  mag     e_KsVISTA    ?=- Ks VISTA magnitude error
  962- 967  F6.3  mag       W1         ?=- W1 magnitude
  969- 973  F5.3  mag     e_W1         ?=- W1 magnitude error
  975- 980  F6.3  mag       W2         ?=- W2 magnitude
  982- 986  F5.3  mag     e_W2         ?=- W2 magnitude error
  988- 993  F6.3  mag       W3         ?=- W3 magnitude
  995- 999  F5.3  mag     e_W3         ?=- W3 magnitude error
 1001-1006  F6.3  mag       W4         ?=- W4 magnitude
 1008-1012  F5.3  mag     e_W4         ?=- W4 magnitude error
--------------------------------------------------------------------------------

Acknowledgements:
    Richard Culpan, rick(at)culpan.de

(End)                                        Patricia Vannier [CDS]  12-Mar-2022
