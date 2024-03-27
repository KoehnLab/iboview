// How to use this:
//  - setup your scene in IboView (orbitals, colors, modes, etc)
//  - open this (js) file in a text editor
//  - adjust the target file name to your liking
//    (see comments at end)
//  - select the text in the entire file and copy it to the clipboard
//    (as text). Typically, Ctrl+A then Ctrl+C. Even Notepad can do it.
//  - switch to IboView and use Edit/Exec Script...
//    in the menu (or just press Ctrl+Shift+V).
//    This will execute the script.
//  - After IboView does the rendering and saving, you should find
//    a bunch of .png files named as you chose in the target directory.



// turn off image cropping (to make sure that all frames have the same size)
view.crop_images = false;

// iterate through all the frames loaded. Indexing starts at 0 and goes
// to num_frames() - 1.
for (var iFrame = 0; iFrame < app.num_frames(); ++iFrame) {
   // switch to frame iFrame.
   app.set_frame(iFrame);

   // make a file name (see below for comments)
   var filename = format("/tmp/claisen_f{0}.png", fmti("%04i",iFrame));

   // save whatever we currently see under the filename created above.
   view.save_png(filename);
}

// Notes on the creation of the file name:
//    var filename = format("/tmp/claisen_f{0}.png", fmti("%04i",iFrame));
//
//    - In windows, you may need to 'escape' backslash characters as
//      directory separators by doubling them. That means instead of writing
//      "x:\data\frame{1}.png", you might have to write
//      "x:\\data\\frame{1}.png". The reason for this is that some
//      character combinations starting with '\' have a non-literal meaning
//      in C-like languages (e.g., '\n' is the escape sequence for newline).
//      Alternatively, use front-slashes as directory separators in Windows,
//      too. It is uncommon, but Windows will happily take a filename
//      like "x:/data/frame{1}.png" to mean "x:\data\frame{1}.png".
//
//    - the 'format' function gets a template (here "/tmp/claisen_f{0}.png")
//      as first argument
//    - and then replaces '{0}' in the template by the
//      first argument after the template (here "fmti("%04i",iFrame)")
//    - ...if present, replaces '{1}' in the template by the
//      second argument after the template
//    - ...if present, replaces '{2}' in the template by the
//      third argument after the template...
//
// We here give only one argument, "fmti("%04i",iFrame)". This
// means 'format an integer (i) in such a way that the total length is 4,
// and filling to the left with zeros if the integer is not large enough'
// (e.g., it would format the integer '2' as string '0002'. In contrast,
//  fmti("%06i",2) would result in the string '000002' and
//  fmti("%4i",2) (without the 0) would result in '   2' (with spaces instead
//  of zeros))
