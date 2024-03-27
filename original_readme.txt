* Updates & Information
   - For updates and information, see http://www.iboview.org

* Prerequisites
  * Graphics
    - OpenGL 4.0-capable graphics adapter
    - Working OpenGL 4.0 drivers
  * Libraries
    - Either QT4 or QT5 development libraries. This IboView uses the QT
      modules: 
        
        core gui widgets opengl script svg
        
      Depending on the system, some of those may need to be installed
      separately. E.g.:
    
      - For QT4 on Fedora/Red Hat, the following works for me:

          dnf install qt4-devel
                  
      - But for QT5 development libraries, qtsvg and qtscript need to
        be installed separately:
      
          dnf install qt5-devel qt5-qtsvg-devel qt5-qtscript-devel 
    
    - On SuSE: A separate installation of the package 'glu-devel' via
      YaST may be required if compilation fails with a message about
      missing 'GL/glu.h'.

  * Mathematical libraries
    + In IboView 2019+, an explicitly provided BLAS/LAPACK is no longer
      strictly required. IboView comes with an internal packaging of the
      C++ Eigen matrix library (https://eigen.tuxfamily.org/index.php?title=Main_Page)
      and a minimal subset of f2c'd LAPACK functions (http://www.netlib.org/lapack/).
      
      This version is normally sufficient if IboView is used for chemical
      analysis only, but for computing wave functions internally via MicroScf,
      it may imply severe performance penalties.
   
    + For better performance and stability, we recommend using 
      the Intel Math Kernel library (MKL) or another working 64bit
      BLAS/LAPACK (e.g., OpenBLAS).
      If the MKLROOT environment variable is set (this may require a
      "load module mkl" first), the main project file should be able
      to pick up and use MKL automatically. Otherwise, one may overwrite
      it manually at the strat of main.pro (some examples are provided
      in there) or set it as environment variable explicitly:
      
      MKLROOT=/opt/intel/composerxe/mkl qmake-qt4 main.pro && make -j 12

* How to compile
   - To compile the program: In the main directory (the directory with the
     file "main.pro") execute:
   
        mkdir build && cd build
        qmake-qt4 ../main.pro && make -j 12 && make install
        cd ..

     Note:
     + "make install" just copies the iboview executable from the build
       directory to the main directory. It does not attempt a system-wide
       installation, and it does not require root/admin privileges.
     + If `qmake-qt4` is not available, try `qmake-qt5`, or just `qmake`.
       If neither work, you may need to locate a suitable installation via:
       
         locate qmake-qt4
         locate qmake-qt5
         locate qmake
      
       If these still do not work, make sure that the required QT development
       libraries are installed (see Prerequisites).
     + If MKL is installed and MKLROOT is correctly set (see Mathematical Libraries),
       this procedure should pick up MKL automatically. If it cannot find MKL,
       it will use an hacky internal BLAS/LAPACK subset, based on C++ Eigen
       and f2c'd LAPACK functions. To prevent that, you can modify the entries
       at the start of main.pro to use other BLAS/LAPACK lines.

   - If the compilation succeeded, try
   
       'cd example-data && ../iboview feco3no_ibos_exp2.js'
      
     to test the IboView program.


* Input data
  * The program can read the following file types:
    * Natively:
      - .xyz files (with one or multiple geometries, with or without
         energy and gradient information)
      - molden files
         (one geometry only)
      - Molpro xml files (.xml) from Molpro 2012 and Molpro 2015
         (one geometry only)
      - Its own script files (.js)
      - Turbomole "coord" and "gradient" files
    * With support by external converters:
      Additionally, IboView can be setup to use several conversion programs
      automatically, so that other types of inputs can be loaded like
      natively supported files:
      - Orca .gbw files
      - Turbomole control files
      The converters must be setup under: Edit/Preferences/External Tools
  * Notes on IRC calculations (or similar series of calculations)
    - We provide a set of examples for running orca or Molpro on series
      of geometry frames (e.g., from an IRC) to create and store wave
      function files. They should be easily adaptable to other programs.
    - These are provided in the `irc-scripts` subfolder of IboView.
    - IboView's built-in tools to load, align, and re-order frames, and
      to export them to a single .xyz file (load data, Press Ctrl+F)
      can be very useful to create the input data for such calculations.
      E.g., for merging TS+both sides of an IRC trace into a single .xyz
      file in consistent alignment and order.
  * Using with Molpro:
    - Do a regular calculation.
    - Molpro 2015: Export last orbital set to xml files via

        {put,xml,filename.xml; nosort; keepspherical}

    - Molpro 2012: Export xml files via

        {put,xml,filename.xml; nosort}

    - Complete example inputs are provided in the example-data directory.
      (see, e.g. trans_b_SN2_intermediate.inp).
  * Using with Turbomole:
    * Simple method (make IboView do conversions automatically)
      - If you have a working Turbomole installed locally (i.e., on
        the computer running IboView), you can just set up the
        path to "tm2molden" in IboView under:
        
        Edit/Preferences/External Tools.
      - Afterwards IboView will run `tm2molden` automatically if you
        attempt to load `control` files (e.g., by dragging and
        dropping them into the IboView main window).
      - IboView can read turbomole coord/gradient files natively
        in any case (i.e., does not actually need Turbomole's
        t2x conversion programs)
      - Details: For the Turbomole control file interface, IboView calls Turbomole's
        file conversion program tm2molden, and this program needs to be set up under
        Edit/Preferences/External Tools. The default is
        "$(TURBODIR)/bin/$(TURBOMOLE_SYSNAME)/tm2molden", using the TURBODIR and
        TURBOMOLE_SYSNAME environment variables. This should work in environments
        in which Turbomole itself is setup to run (e.g., on my system, I'd run
        export TURBODIR=/opt/TURBOMOLE/ && source $TURBODIR/Config_turbo_env
        where $TURBODIR points to my local Turbomole installation. In your system
        this may be different.). 
        If IboView is started in this environment afterwards, it should be able to
        load control files directly.
    * Complicated method (no access to turbomole itself):
      - Run your calculation, and export its results to molden format
         via Turbomole's `tm2molden` program.
      - Read the Molden files with iboview.
      - You can also use `t2x` to just export geometry data to .xyz file
         format (for example, to track geometry optimization progress)
  * Using with Orca:
    * Simple method (make IboView do conversions automatically)
      - If you have a working orca installed locally (i.e., on
        the computer running IboView), you can just set up the
        path to `orca_2mkl` in IboView under:
        
        Edit/Preferences/External Tools.
      - Afterwards IboView will run `orca2_mkl` automatically if you
        attempt to load `gbw` files (e.g., by dragging and
        dropping them into the IboView main window). Note that for
        large files this can take a while, during which IboView will be
        unresponsive.
      - Details: For the Orca .gbw file interface, IboView calls Orca's orca_2mkl program.
        For this to work, the path to the program has to be set up in Edit/Preferences/External Tools.
        E.g., on my system I enter "/opt/prog/orca_4_0_1_2_linux_x86-64_shared_openmpi202/orca_2mkl",
        where /opt/prog is where I installed Orca to. This should work on Windows, too
        (by entering the corresponding orca_2mkl.exe file path).
        + Note that for this to work on linux, orca itself needs to be correctly setup so that it
          can run in the same environment---for the MPI version of Orca that means that a compatible
          OpenMPI is setup in the environment variables.
        + @Linux: In case you run Orca on a cluster and cannot get the MPI version running locally on
          the machine on which IboView runs, the easiest way to get this working might be to install a
          serial Orca (non-MPI) additionally and use the orca_2mkl of that one (note that gbw files of
          different versions might not be compatible, so if one is changed, the other one should be
          changed, too).
    * Complicated method (no access to orca itself):
      - Use "orca_2mkl <basename> -molden" to export your computation
        result to molden format (e.g., if you have a wave function file
        called "wheee.gbw", you would call "orca_2mkl wheee -molden". Note
        that the .gbw file extension is not given to orca_2mkl.
        `orca_2mkl` then generates a file called wheee.molden.input)
      - Read the Molden files with iboview.
  * Using with Molcas:
    - Import of .molden files generated by Molcas should now work.
  * The program can read multiple input files at once. If invoked via

         iboview file000.xml file001.xml file002.xml ...

    (or iboview dir/*.xml) it will load the given files and attempt
    to sort the contained orbitals by overlap, possibly flipping their
    phase as required. Properties (e.g., color) of associated orbitals
    will also be linked.
    The different files can then viewed by changing the frame number
    in the "Data Sets" page. Pressing Ctrl+T will render iso surfaces
    for the selected orbitals in all frames at once (instead of on
    demand).
  * Other input methods:
    - The program accepts all supported types of files via drag & drop
      (e.g., drag files from your file browser into an open IboView
      window to open them)
    - The state of the visualization can be copied to clipboard via
      Ctrl+Shift+C, and later restored via Ctrl+Shift+V (Hotkeys for
      Edit/Copy State and Edit/Exec Script). The states take the form of
      .js scripts and can be saved, edited, or combined externally in a
      text editor if desired. For example, it is possible to only
      restore the viewing angle or orbital colors by removing all
      unrelated commands from the script before executing it via
      Ctrl+Shift+V.
    - .xyz Files can also be pasted from Clipboard and need not be
      saved to disk: If the text of an .xyz file is in the clipboard,
      Pressing Ctrl+Shift+V in IboView (Edit/Exec Script) will
      load it.

  * Notes on security:
    - Do not load data files from people you do not trust. The program has
      not undergone ANY security-related code-review and most likely
      contains many bugs which can be exploited to execute arbitrary code.

    - In particular, you should UNDER NO CIRCUMSTANCES execute state
      scripts (via loading .js/.chai files or using "Exec Script
      (Clipboard)") from people you do not trust or which you have
      not yourself determined to be safe.
      State scripts are PROGRAMS, and could do dangerous things. The
      authors take no responsibility for what happens when they are
      executed, and in particular DO NOT guarantee that the scripts are
      sufficiently sandboxed to prevent them from affecting your system
      or data. You have been warned.

* Using and controlling the program
  * Controls: Data sets
    - Double-click on an orbital in the "Data Sets"-page or
      "Find Orbitals"-dialog to toggle rendering it
    - Single click on a /visible/ orbital in either the "Data Sets" page
      or the 3d-view in order to change its properties in the "Render:
      Orbital" page.
    - Right-clicking on an object (orbital, atom, or bond) in the
      3d-view brings up a context menu for this object.

  * Controls: Finding orbitals
    - For getting an overview of the orbitals, it may be helpful to
      reduce the iso resolution to 10 pt/A. This speeds up surface
      construction by approximately a factor of 5...10.

    - To find orbitals localized on specific atoms, select one or more
      atoms (by clicking them and holding Shift), then right click on one
      of the selected atoms, and select "Find Orbitals..." in the context
      menu.

    - In reaction mechanisms: In DataSets/Frames, you can get an idea of
      which orbitals *change* along the currently loaded frames by
      turning the "Track Orbital" dial. The white dotted line in the
      curve view will then adjust to show how the current orbital number
      changes along the reaction path (it shows the root mean square
      change of the orbital's atomic charges compared with the first
      frame).

      If the line stays flat, it means that the orbital's charge does not
      move between atoms. If it goes up along the reaction path, it means
      that the orbital is distorted during the reaction, possibly
      breaking or forming new bonds.

      If you find an orbital which changes, press the "Show #.." button
      to toggle rendering it, and press Ctrl+T, to render its iso-
      surfaces for all frames. Afterwards, you can use the "current
      frame" dial to see how the orbital transforms during the reaction.
    
    - A summary table of which orbitals change most along the loaded frames can
      also be accessed under "View/Find Actively Reacting Orbitals...". This
      provides a list of the orbitals with the largest atomic partial charge
      displacements along the loaded frames. Double clicking inside the table
      toggles viewing an orbital (like in the other "Find..." tables).

  * Controls: Moving in the 3D view
    - Hold right mouse button and move mouse: rotate camera
      ...do so while holding the Shift key: rotate around current
          position instead of coordinate origin.
    - Hold right mouse button and turn wheel: zoom in/out
    - Hold right+left mouse button and move: translate (move picture)
    - Hold right+left mouse button and turn wheel: rotate view around mouse
    - Just turn wheel: Currently do nothing.
    
    - Move between frames:
      + Use frame controls on Main Window / Data Set/ Frames
      + Or use keyboard:
        - First click on the 3D view to give it focus (e.g., click
          on an empty space)
        - Then use:
          Shift + <Left Arrow>           : switch to previous frame
          shift + <Right Arrow>          : switch to next frame
          Shift + Ctrl + <Left Arrow>    : switch to frame several indices back
          shift + Ctrl + <Right Arrow>   : switch to frame several indices further
      + Or use mouse wheel:
        - First click on the 3D view to give it focus (e.g., click
          on an empty space)
        - Then use:
          Shift + <Mouse Wheel>    : switch between frames

  * Controls: Selecting and grouping atoms in 3D view
    - Atoms can be selected by single-clicking them. Multiple atoms
      can be selected by holding Shift while selecting them (this
      toggles their selection status).
      
    - Atoms within a rectangle can be selected by opening a selection rectangle
      (click on one part of the molecule, hold left mouse button, drag somewhere
      else, release). This currently yields no visual feedback (not yet
      programmed...), but selects all the topmost visible atoms in the given
      rectangle.

    - Selection groups:
        + Selections can be saved by pressing Ctrl+1, Ctrl+2, ... to create
          selection groups
        + Selection groups can be restored by pressing 1, 2, ....
        + To clear a selection group, make an empty selection (e.g.,
          by single clicking on an empty space next to the molecule), and
          press Ctrl+#n (with n=0,1,...) to assign the empty selection to
          the group you wish to clear.
      Some analysis programs (e.g., the EOS oxidation state analysis) use
      these selection groups to define the desired fragmentation of the
      system for analysis.

  * Rendering Orbitals & Iso-Surfaces
    - If you see artifacts in places where many orbitals intersect,
      increase the number under "Depth Layers" in Render: Orbitals/Surfaces.
      This specifies the maximum depth complexity of transparent surfaces
      in the program (increasing this number will slow down interactive
      rendering, but not iso-surface tracing).
    - To speed up the construction of iso-surfaces, decrease the number
      in Render: Orbital / Surfaces / Resolution (1/A). Iso-surface time
      scales with r^3. The default setting of 20/A is recommended for
      publication quality graphics (although other programs tend to use
      much smaller resolutions), for survey studies 10/A is sufficient.

      
* Analysis methods implemented in IboView
  * IAOs/IBOs
    - TODO: describe 2013 IAO paper, 2015 Curly arrow paper, 2019 Ni paper for
      virtual valence IBOs, and 2021 relativistic IBO paper (for fragment-localization
      and writeup of 2014' update of IAO construction)
  * Bonds orders
    (not published yet)
  * Effective Oxidation State (EOS) Analysis
    * TODO: c/p or describe EOS analysis core idea
    * Notes on EOS analysis outputs:
      - When the EOS analysis gives a reliability of the assignment of "50%", that really means that there is no reliability at all (Salvadore defined the r-number like this). If something like this happens, what it means is that the system does not have a real oxidation state, because it features very strong covalent bonding, and the ionic picture breaks down (at least according to the input wave function).

      - However, that itself is often a matter of interest to experimentalists in some cases. So this may not be necessarily bad, especially if unexpected. In this case I would recommend having a look at the orbitals themselves to clarify the bonding situation and to describe the actual bonding motive as a new discovery. At least this has worked for me and Johannes a few times...

      - Note that the EOS assignment reliability is almost entirely influenced by the input wave function---in this the orbitals are actually not localized in the MO localization sense in that where one orbital goes, influences the position of the other. The reliability rating cannot be transformed away. So multiple very close EOS assignments is normally something I would regard as indicator of an actual physical phenomenon.



* Scripts
  - The program can be controlled by scripts of Javascript format, but
    the documentation is not written yet.
  - Some example scripts are provided in the "example-scripts" subdirectory
    provided with IboView.


* Known problems
  - On Mac, the alpha channel has to be disabled completely because
    MacOS uses it for blending the view itself (?!). Run iboview with
    -disable-alpha and -disable-backfaces commandline arguments for
    workarounds.

  - The embedded MicroScf program is still quite limited (open-shell
    not very clever, no ECPs slow and no gradients, optimally cached/
    semi-direct SCF, no convergence stabilizers beyond level shifts and
    DIIS, only a few xc functionals, ...).
    Additionally, the front-end in IboView does not allow specifying
    several advanced options, and the memory-for-method computation
    is not properly taking account all factors.
    Will be fixed later.
    The default DFJX-RKS/PBE should be good for tracing organic compounds,
    however.

  - Also, internally computed wave functions can currently not be stored
    and reloaded.

* Meta    
  * Other used software libraries
    The following software is distributed in source form together with
    iboview:
    - glew: The OpenGL Extension Wrangler Library, libGLEW.
      iboview employs GLEW in order to access OpenGL functions.
      See http://glew.sourceforge.net/ and
      http://glew.sourceforge.net/credits.html for credits.
      GLEW is licensed under the Modified BSD License, the Mesa 3-D License
      (MIT License), and the Khronos License (MIT License).

    - pugixml:
      iboview employs the pugixml library (http://pugixml.org),
      pugixml is Copyright (C) 2006-2014 Arseny Kapoulkine.
      pugixml is distributed under a MIT license. See pugixml/pugixml_readme.txt

    - QPropertyModel.h/.cpp:
      A class for easily turning any QObject-derived subclass with
      properties into a one-row model (cgk: i.e., a QAbstractItemModel/
      QDataWidgetMapper interface for QT properties, which allows
      easy syncing of UI controls with data. It's pretty neat.
      Look it up.).
      Copyright 2013 - Harvey Chapman <hchapman@3gfp.com>
      See https://gist.github.com/sr105/7955969

      QPropertyModel is licensed under the Creative Commons Attribution-
      ShareAlike 4.0 International License. To view a copy of this
      license, visit
      http://creativecommons.org/licenses/by-sa/4.0/deed.en_US.

    - shader/pixel5_fakeAA2.glsl:
      This is a FXAA ("fast approximate anti-aliasing") filter. The
      version here is a slightly modified version (added alpha channels)
      of the GLSL version ported by Sebastian F. Mazza
      (http://code.google.com/p/blubbengine2/) from Timothy Lottes'
      (http://timothylottes.blogspot.de/) original HLSL version. My
      understanding is that the original HLSL FXAA code was made public
      domain by nvidia (an earlier blog post of Lottes said so), the port
      of Sebastian F. Mazza seems to be released under "other open source
      license" according to http://code.google.com/p/blubbengine2/ (the
      actual license was not included).

    - A distance-field pixelized version of the Liberation Sans font is
      embedded into the program. Liberation Sans is released under the
      SIL Open Font License.
      See https://fedorahosted.org/liberation-fonts/

    - fn_LiberationSans.h and shader/pixel_text.glsl:
      The former is an interface file for the Liberation Sans
      font as generated with a modified version of makefont of
      Freetype GL (http://code.google.com/p/freetype-gl/).
      The pixel shader pixel_text.glsl is a modified version
      of the shader distance-field.frag employed in freetype-gl's
      demo-distance-field.c.

    - vector_math.h:
      A modified version of vector_math.h, Copyright (c) 2007,
      Markus Trenkwalden, from http://www.trenki.net/files/vector_math.h

  * On Forking and Contributing
    - Do not fork the program. Doing so will invoke my wrath.

      You can freely adjust the program to your own/own group's/own
      institute's needs, but please do not distribute the modified
      program. Especially not via freely accessible web pages.

    - I am ready to accept contributions to the main code base.
      Particularly helpful would be contributions which improve the user
      experience, e.g., support for additional import formats.

      However, there are restrictions to code contributions:
      - Changes will only be incorporated into the main code base if you
        grant me, Gerald Knizia, an unrestricted and unconditional license
        to use, modify, distribute, and re-license the contributed code.
        (in short: if the contributed code will create potential legal
        problems in the future I do not want it).
      - If you make changes which you wish to incorporate into
        the main code base, these /#have to be coordinated with me#/.
      - Changes are unlikely to be accepted if they involve serious
        maintanance- or deployment burdens (e.g., by introducing 3rd party
        libraries with license problems or which are painful to build on
        MSVC), or involve methods I do not like (e.g., AIM).
      - Changes which are developed or published without coordination
        with me will /*DEFINITELY*/ not be accepted[1]

      [1] The development and maintanance of IboView was and is highly
      non-trivial. Please consider that whatever you are currently doing
      might also have been sitting on my todo list long before I ever
      heard of your publication. However, my resources are limited; I
      am thus reserved about possibily promoting other people's work via
      IboView which I could not do myself solely because I spent the
      resources on IboView itself. For this reason developments should be
      coordinated with me, and some of them might be acceptable only if
      done in cooperation.

  * Change History:
    * v20210521-RevA (2022-01-23):
      Minor bugfix update (via patches to the v20210521 code base):
      - Add the forgotten .rc definition which caued IboView's icon to be missing
        from the exe file. Clearly the most important update :)
      - Fix an error in tabulated free-atom configurations, which prevented
        MicroScf from generating SCF initial guesses for Pd (palladium) atoms
        (thanks to Prof Russel Hughes for bringing this to my attention).
        The error did not affect any other elements than Pd or any other
        tasks than generating SCF initial guesses---in particular, it did not
        affect IAO/IBO analysis of wave functions imported from other programs.
      - Work around a potential issue with MSVC's OpenMP implementation which
        could get it to hang in DFT grid instanciation of large molecules with
        many threads.
      - Clean up qmake configuration (main.pro)
    * v20210521:
      Feature updates:
      - Various bugs/errors in IboView have been fixed. In particular, an annoying previous crash when loading molden files and trying to geometry-align them. Ctrl+F Ctrl+F on loaded files should be safe now.
      - IboView should deal more gracefully with very long IRCs. Additionally, features to compress very long IRCs have been added (see Ctrl+F/Delete Frames Closer Than ...).
      - IboView can now reload files (F5). This closes the current file and reopens it (e.g., for following running geometry optimizations).
      - IboView can now compute bond orders. The methodology paper is not published yet, though.
      - IboView can now compute iso-density surfaces, and spin-density surfaces (Main Menu/Compute/Density...).
      - IboView can now compute Effective Oxidation States (Main Menu/Compute/Oxidation States... see below).
      Major changes to the 2015 version:
      - IboView will now remember settings, positions/sizes of windows and dialogs, and locations of last opened and saved files.
      - Dragging & dropping a file into IboView, or using File/Import Files... will no longer close the current file automatically. The previous one has to be explicitly closed (File/Close or Ctrl+W)
      - The development version will no longer load multiple frames with different numbers or types of atoms. Old files need to be opened before opening new ones with different molecules. Plan to fix this before 2019 release, but at this moment the program will most likely crash if opening is attempted.
      - If wave function information is available, IboView will now render bond topology based on computed bond orders (including rendering partial bonds and multiple bonds), rather than infer bonds based on atomic distance only. The previous geometry-based bonds can be re-established by using Edit/Reset Bond Lines (Ctrl+Alt+R), or by changing the Main Window/Render: Geometry/Bonds/Source setting.
    * v20150427:
      - Fix a problem with Windows vs Unix line endings, which prevented
        IboView from reading Unix-generated .molden-Files in the Windows
        version.
      - Changed default IBO exponent from 4 to 2.

    * v20150424:
      - .molden files from Orca and Molcas are now supported.
        A sanity check was added to IBBA, to diagnose broken imports
        (via non-orthogonal occupied orbitals)
      - Fixed a problem with .molden files from Turbomole, which could
        lead to an application freeze.
      - Settings regarding window sizes and file histories are now remembered
        across application starts.
      - Added support for automatic execution of a specified script file
        on startup.
      - Properties of elements and atoms (e.g., colors, sizes, covalent radii,
        etc) can now be changed via scripts.

    * v20150214:
      - Initial public release



#+ kate: indent-width 2
