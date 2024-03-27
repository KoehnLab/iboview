// A simple example on how to change atom colors for highlighting purposes.
// How to use: c/p script text and Ctrl+Shift+V in IboView to execute
// (or drag and drop .js file into IboView).
//
// For additional details: see 
// - change_atom_or_element_properties_00_reference_info.js
// - change_atom_or_element_properties_01_combining_with_loops_and_arrays.js

var atoms_to_highlight = [2,3,4,5,6];
// ^-- to see which ones: use Ctrl+N to toggle atom numbers.

for (var i = 0; i != atoms_to_highlight.length; ++i) {
   var iAt = atoms_to_highlight[i];
   doc.atom_options(iAt).color = 0x903030;
   // ^- this is the color to change, in HTML RRGGBB format (hexadecimal from 00 to ff)
}
doc.update_views();


// if you want to change other aspects of the atoms:
// apart from '.color', there are:
//
// bond_color
// draw_radius
// vdw_radius
// draw_name  (defaults to element name)


// note:
// doc.reset_atom_options(2);
// resets all changes to defaults for atom 2
