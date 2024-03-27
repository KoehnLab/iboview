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


var mos = [];

mos.push([8, 0x99ff66a6, 0x99ffa666]); // MO number, color for positive side, color for negative side
mos.push([22,0x99ff66a6, 0x99ffa666]); // (these can be obtained by setting up all orbitals, and hitting Edit/Copy State)
mos.push([23,0x9966a6ff, 0x99a666ff]);

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


