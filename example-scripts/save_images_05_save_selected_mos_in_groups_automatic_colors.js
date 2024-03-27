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


var mo_groups = []; // will accumulate the groups of MOs to display together (including colors).
var icolor_phase = 0; // will be incremented after each MO is added

function add_mo_group(iMos) {
   // given a list of MO indices, add a visualizaiton group for them. Assign colors to each MO.
   var mos = [];
   for (var iiMo = 0; iiMo < iMos.length; ++ iiMo) {
      var h = 45.0 * icolor_phase;
      var a = 0.6;
      var v = 1.0;
      var color_plus = hsva(h, 0.7, v, a);
      var color_minus = hsva(h, 0.3, v, a);

      mos.push([iMos[iiMo], color_plus, color_minus]);
      icolor_phase += 1;
   }
   mo_groups.push(mos);
}


// add MOs 8, 22, 23 to a single visualization group (shown together); Take assign colors automatically.
add_mo_group([8,22,23]);
add_mo_group([19,20,14]);

for (var iMoGroup = 0; iMoGroup < mo_groups.length; ++ iMoGroup) {
   // get a reference to the MO lists for the current group.
   var mos = mo_groups[iMoGroup];
   for (var iFrame = 0; iFrame < doc.num_frames(); ++iFrame) {
      app.set_frame(iFrame);
      var frame = app.frame();

      // switch on the MOs in the current group.
      for (var iiMo = 0; iiMo < mos.length; ++ iiMo) {
         var mo = mos[iiMo];
         app.show_mo(mo[0], mo[1], mo[2]);
      }

      var filename = format("/tmp/frames/claisen_group{0}_f{1}.png", fmti("%02i",iMoGroup), fmti("%03i",iFrame));
      print("saving: " + filename);
      view.save_png(filename);

      // switch them off again.the MOs in the current group.
      for (var iiMo = 0; iiMo < mos.length; ++ iiMo) {
         var mo = mos[iiMo];
         app.hide_mo(mo[0]);
      }
   }
}

// kate: syntax javascript;


