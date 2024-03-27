view.set_size(933,811);
view.camera.set_pos(-4.58031,-16.10371,10.94035);
view.camera.set_dir(0.22902,0.80519,-0.54702);
view.camera.set_vup(0.11581,0.53542,0.83661);
view.camera.set_zoom(1.00000);

// when making videos, you probably want to turn of auto-crop (as otherwise all
// the exported frames will have different sizes) and saving alpha channels
// (most video software cannot cope with the transparent backgrounds).
// This can be done by pressing the corresponding buttons on the UI, or including
// this in the scripts:
view.save_alpha = false;
view.crop_images = false;


var mos = []; // will accumulate the MOs to display (including colors).
var icolor_phase = 0; // will be incremented after each MO is added

function add_mo(iMo) {
   // make a color for the positive and negative sides. Colors are 32bit unsigned integers of
   // the following format:
   //     0xAABBGGRR,
   // where AA is the hexadecimal opacity betwee 0 and 0xff (=255), and the RR (=red) GG (=green)
   // and BB (=blue) formats work similarly. E.g. 0xff0000ff is a 100% red at full opacity, and
   // 0x7f0000ff is 100% red at half opacity.
   //
   // IboView includes the following functions for color manipulation:
   //   - hsva(h,s,v,a): Make a RGBA color based on the hue/saturation/value color model.
   //     h is a phase angle between 0 and 360., s, v, and a are floats in interval [0.,1.].
   //     'value' is a form of expressing color brightness. 1.0 means full brightness of the
   //     given color (unlike in HSL, this does not go to white).
   //     In IboView, default values of saturation are 0.6, value 1.0, and opacity (=alpha) 0.6 (60%).
   //   - irgb: (invert RGB). This exchanges the red and blue channels of a color, to account for
   //     differences in HTML and OpenGL color models.

   var h = 45.0 * icolor_phase;
   var a = 0.6;
   var v = 1.0;
   var color_plus = hsva(h, 0.7, v, a);
   var color_minus = hsva(h, 0.3, v, a);

   mos.push([iMo, color_plus, color_minus]);
   icolor_phase += 1;
}


// add MOs 8, 22, 23; Take assign colors automatically.
var imos = [8,22,23];
for (var iiMo = 0; iiMo < imos.length; ++ iiMo)
   add_mo(imos[iiMo]);


print(format("mos.size: {0}", mos.length));
for (var iiMo = 0; iiMo < mos.length; ++ iiMo) {
   var mo = mos[iiMo];
   var iMo = mo[0];
   app.show_mo(mo[0], mo[1], mo[2]);

   for (var iFrame = 0; iFrame < doc.num_frames(); ++iFrame) {
      app.set_frame(iFrame);
      var frame = app.frame();
      var filename = format("/tmp/frames/claisen_mo{0}_f{1}.png", fmti("%02i",iMo), fmti("%03i",iFrame));
      print("saving: " + filename);
      view.save_png(filename);
   }
   app.hide_mo(mo[0]);
}

// kate: syntax javascript;


