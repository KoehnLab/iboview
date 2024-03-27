app.load_file("trans_b_sn2_intermediate.xml");
// add partial bond lined for bonds being made/broken
doc.add_bond(1,3,"dotted|gray"); // C-C banana bond
doc.add_bond(1,5,"dotted|gray"); // C-O bond

view.camera.set_pos(-7.64953,18.47107,-0.79091);
view.camera.set_dir(0.35682,-0.93360,0.03287);
view.camera.set_vup(-0.90487,-0.35416,-0.23617);
view.camera.set_zoom(0.82270);

doc.show_mo(75,0x99ffa666,0x99ff66a6);
doc.show_mo(79,0x9966ffa6,0x99a6ff66);
doc.show_mo(81,0x99a666ff,0x9966a6ff);

// add some more orbitals. we here set the colors manually
// using the hsva(h,s,v,a) function, which assembles a color
// from hue (0...360), saturation (0..1), value (0..1), and alpha(0..1)
// alpha specifies opacity. 0.6 corresponds to the 0x99000000 used above.
var h = 0.0;
var sm = 0.3;
var sh = 1.0;
var a = 0.6;
// h = 155.0; doc.show_mo(75,hsva(h,sh,1.0,a),hsva(h,sm,1.0,a));
// h =  90.0; doc.show_mo(79,hsva(h,sh,1.0,a),hsva(h,sm,1.0,a));
// h =  45.0; doc.show_mo(81,hsva(h,sh,1.0,a),hsva(h,sm,1.0,a));
h =   0.0; doc.show_mo(83,hsva(h,sh,1.0,a),hsva(h,sm,1.0,a));
view.save_png("trans_b_sn2_intermediate.png");

// kate: syntax javascript;
