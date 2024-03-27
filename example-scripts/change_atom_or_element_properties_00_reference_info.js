// coding: utf-8

// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
// File Description:
// ─────────────────────────────────────────────────────────────────────────────
// - This script demonstrates how to change the properties (e.g., drawing size,
//   atomic/element labels, colors) of elements and/or individual atoms.
//
// - By default, all of these settings apply only to the currently running
//   IboView instance. That is, they will be reset to their defaults if IboView
//   is restarted. So you can experiment with these as desired without breaking
//   the defaults.
//
// - The settings can be either added to js script files for specific data
//   files, or executed by copying the script text (Ctrl+A, Ctrl+C) and pasting
//   it in IboView via Ctrl+Shift+V (that is Edit/Exec Script...).
//
// - If you would like to make any of these settings permanent (i.e., set them
//   up as new global defaults), you can create a new .js IboView script file
//   containing the settings you would like to change, and tell IboView to
//   load and execute this file automatically on startup under
//   Edit/Preferences/General Settings
// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

// ░░░░░░ global scales ░░░░░░

// you can change the scale of *all* atoms/bonds/labels by adjusting the
// following view scales, which are given in percent. Changing these by script
// is equivalent to turning the dial controls under "Render: Geometry"
// (atoms/bonds) or spin control (labels) in the IboView main window
view.atom_scale = 100.0;  // scale of all atom balls; in %, relative to IboView default.
view.bond_scale = 100.0;  // scale of all bond lines; in %, relative to IboView default-times-atom_scale.
view.label_size = 100.0;  // scale of atom labels for elements (Ctrl+E) and/or atom numbers (Ctrl+N)


// ░░░░░░ element properties ░░░░░░

// change colors of elements hydrogen (all Hs), carbons (all Cs), and oxygens (all Os).
doc.element_options(1).color = 0xff0000;  // make hydrogens red (element 1)
doc.element_options(6).color = 0x00ff00;  // make carbons green (element 6)
doc.element_options(8).color = 0x0000ff;  // make oxygens blue (element 8)
// ^-- colors are given as a composite hexadecimal number (prefix '0x') in
// format 0xRRGGBB, where RR/GG/BB are two-digit hexadecimal numbers (00,...,ff)
// controlling the red/green/blue channels

// you can also change the element label from their default element symbols.
// E.g., the next line changes the element labels of carbons from "C" to
// "carbon" (yes, this will look silly)
doc.element_options(6).draw_name = "carbon";


// ░░░░░░ atom properties ░░░░░░

// Next to changing the properties of *elements* in general, you can also
// change properties of individual *atoms*, which are identified by number
// (use Ctrl+N / Ctrl+E to toggle atom numbers/element labels)

// E.g., change color and label of *only* atom 3, regardless of which element
// it is. This will not affect other atoms, even if they share the same element.
doc.atom_options(3).color = 0x303030;
doc.atom_options(3).draw_name = "Nu";


// ░░░░░░ re-rendering the scene with updated settings ░░░░░░

// after changing the properties, an explicit update of the scene is required
// before they become visible. update_views() forces this (otherwise, moving in
// the scene or resizing the window will also do it).
doc.update_views();


// ░░░░░░ restoring default atom or element properties ░░░░░░

if (false) {  // <- a block in "if (false)" is not executed. Change to "if (true)" to see what the code does.
   // use these to reset the atomic options or element options to their defaults
   // (undoing the changes from above)
   doc.reset_atom_options(3);  // reset atom 3
   doc.reset_element_options(1); // reset element 1 (hydrogen)
   
   // you can also do this in a loop, to reset all elements/atoms
   for (var iElement = 1; iElement <= 100; ++ iElement) {
      print(format("resetting options of element {0} {1}", iElement, doc.element_options(iElement).draw_name));
      doc.reset_element_options(iElement);
   }
   // note: if resetting just element properties, the indvidual atom
   // property overrides will stay intact.
   var frame = app.frame();
   for (var iAt = 1; iAt <= frame.natoms; ++ iAt) {
      print(format("resetting options of atom {0} {1}", iAt, doc.atom_options(iAt).draw_name));
      doc.reset_atom_options(iAt);
   }
   
   // as above, the scene needs to be re-rendered for the changes to become
   // visible.
   doc.update_views();
}




// ▓▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒▒
// ▓
// ▓   ADDITIONAL INFORMATION FOR SPECIALTY APPLICATIONS
// ▓


// ░░░░░░ available element/atom properties ░░░░░░
//
// The full list of properties which can be changed for either individual atoms
// or for all atoms of a given element (see above for examples) is currently as
// follows (see make_properties.py):
//
//    color               # defaults to element default
//
//    bond_color          # defaults to `color`
//
//    draw_radius         # reference radius (comparable to covalent radius); this gets multiplied
//                        # with data derived from view.atom_scale (in %) to control how large atom
//                        # balls (and bonds) are drawn.
//
//    draw_name           # string; specifies the text rendered as "element label" (Ctrl+E)
//                        # for the atoms/element. Defaults to IUPAC chemical element label in default
//                        # case (e.g., "Ar" for argon)
//
//    bond_radius_factor  # decides distance at which bonds are drawn to other atoms
//                        # in geometric bond mode; given as multiple of the covalent radius.
//                        # Defaults to 1.3, meaning that a bond is drawn between A and B,
//                        # if distance(A,B) < 1.3 * rcov(A) + 1.3 * rcov(B).
//
//    vdw_radius          # override for van der Waals radius, given in Angstrom. Default value is currently from
//                        # mixed tables, but will eventually be merged to data from
//                        # 10.1002/chem.201602949 (atom radii from 0.001 e^-/abohr^3 iso surfaces of free atoms).
//    
//    covalent_radius     # override for covalent radius, given in Angstrom. Defaults to a mixture of data derived from
//                        # http://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_%28data_page%29
//                        # and from UFF bond radii



// ░░░░░░ default values for element colors ░░░░░░
if (false) {

   // Note: The default element color scheme of IboView v2015-v2020 is adapted from
   // the "rasmol CPKnew colors" according to
   // http://jmol.sourceforge.net/jscolors/, with some additional substitutions for
   // IboView (e.g., carbon is replaced by 0x999999 (0.6/0.6/0.6)).
   // The concrete default color values are as follows:
   var DefaultElementColors = [0x404040,0xffffff,0xffc0cb,0xb22121,0xff1493,0x00ff00,0x999999,
      0x87cee6,0xff0000,0xdaa520,0xff1493,0x0000ff,0x228b22,0x696969,0xdaa520,
      0xffaa00,0xffff00,0x00ff00,0xff1493,0xff1493,0x696969,0xff1493,0x696969,
      0xff1493,0x696969,0x696969,0xffaa00,0xff1493,0x802828,0x802828,0x802828,
      0xff1493,0xff1493,0xff1493,0xff1493,0x802828,0xff1493,0xff1493,0xff1493,
      0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,
      0x696969,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xa020f0,0xff1493,
      0xff1493,0xffaa00,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,
      0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,
      0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,
      0xdaa520,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,
      0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,
      0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,
      0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493,0xff1493];

// ░░░░░░ changing properties in loops ░░░░░░

   // if you want to, you can also change colors programatically with loops.
   // E.g., like this:
   for (var i = 1; i < 100; ++ i) {
      // set all elements colors to the tabulated values from DefaultElementColors,
      // but with red and green channels inverted ("irgb" changes 0xAARRGGBB to 0xAABBGGRR)
      doc.element_options(i).color = irgb(DefaultElementColors[i]);
   }
   
   // or like this:
   for (var i = 1; i < 10; ++ i) {
      // make rainbow colored elements; hsva = hue (0..360), saturation, value, alpha
      doc.element_options(i).color = hsva(30*i, 1.0, 1.0, 1.0);
   }
   
   // a more complex example of this is given in loop_over_atoms_to_change_properties.js
   
   doc.update_views();
}

