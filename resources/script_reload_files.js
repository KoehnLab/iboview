var open_files = app.loaded_files;
var current_frame = app.get_frame();
app.close_files();
app.load_files(open_files);

function add_quotes(s) {
   return format("\"{0}\"", s);
}

var s = app.loaded_files.map(function(o){return "'{}'".format(o);}).join(", ");
// switch to last frame.. since this would mostly be used to trace geometry optimizations...
// app.set_frame(app.num_frames() - 1);
app.set_frame(current_frame);
app.notify("Reloaded files: {}".format(s));
