// If you copy this into the clipboard (e.g., load into notepad, then do Ctrl+A,
// Ctrl+C) then go to IboView and hit Ctrl+Shift+V (that's paste/execute script
// under Edit) the program will go through all the frames and save each
// picture to a separate .png file in the current directory


// when making videos, you probably want to turn of auto-crop (as otherwise all
// the exported frames will have different sizes) and saving alpha channels
// (most video software cannot cope with the transparent backgrounds).
// This can be done by pressing the corresponding buttons on the UI, or including
// this in the scripts:
view.save_alpha = false;
view.crop_images = false;

for (var iFrame = 0; iFrame < doc.num_frames(); ++iFrame) {
   app.set_frame(iFrame);
   var frame = app.frame();
   var filename = format("ibo-frame-{0}.png", fmti("%04i",iFrame));
   print("saving: " + filename);
   view.save_png(filename);
}
