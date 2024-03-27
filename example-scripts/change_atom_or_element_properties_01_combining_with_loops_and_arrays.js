// see also: change_atom_or_element_properties_00_reference_info.js for a basic
// description of how to change element properties.

// this changes *ALL* oxygens (8 = element number of oxygen; N would be 7, C would be 6, ..., H would be 1, Fe would be 26, etc)
doc.element_options(8).color = irgb(hsva(180., 1.0, 0.8, 1.0));

// ^- note: color has here been generated in HSVA system:
// Hue (0-360), Saturation (0-1), Value (0-1), and Alpha
// (ignored here). The 'irgb' switches the color components
// between OpenGL colors (BGRA), which hsva makes, and
// web colors (RGBA), which atom_options.color takes.
//
// Of course a RGB color (see below) will do the trick, too.


var frame = app.frame();
for (var iAt = 1; iAt <= frame.natoms; ++ iAt) {
   print("now setting color of atom", doc.atom_options(iAt).draw_name, iAt);

   var atom_options = doc.atom_options(iAt);
   if (atom_options.draw_name == "O") {
      // note: colors given in RGB format, #RRGGBB, with each
      // RGB component specified with a hexadecimal number from 00 (black)
      // to ff (full intensity).
      atom_options.color = 0x5050ff;
   }
   
   // of course, if you have many of them, you could just write them in one line:
   if (doc.atom_options(iAt).draw_name == "O") { doc.atom_options(iAt).color = 0x5050ff; }
   // ^- this line does the ones above
   
}
app.update_views();

// note on the loop above:
//  - for compatibility with other programs, the atoms are labelled
//    starting from 1, not from 0.
//  - So in a 12-atom molecule (benzene, for example), the valid
//    atom indices will be 1,2,3,...,10,11,12.


// This can be combined with generic JavaScript (or ECMA-script) features to set properties only of a subset of atoms:

// select some set of atoms to change
var atom_list = [2,5,8,12,37];

// change all atoms from the given list to name 'X',
// and color #40b030, and bigger size, but *ONLY* if they
// are neither Hydrogens nor Carbons.
for (var ii = 0; ii != atom_list.length; ++ ii) {
   var iAt = atom_list[ii];
   var opt = doc.atom_options(iAt);
   print("now setting color of atom", opt.draw_name, iAt);
   if (opt.draw_name != "C" && opt.draw_name != "H") {
      opt.color = 0x60c020;
      opt.draw_name = "X";
      opt.draw_radius = opt.draw_radius * 1.3;
   }
}
app.update_views();


