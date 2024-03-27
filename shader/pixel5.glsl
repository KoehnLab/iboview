layout(location=0) out vec4 color;

// layout(early_fragment_tests) in;
// ^- force depth test enable.
//    otherwise we always accumulate stuff, even if transparent
//    fragment is behind non-transparent fragment.


vec4 calc_base_color(bool FlipSides);

void main(void)
{
   color = calc_base_color(false);

   // set alpha to 1.0 -- the shader currently only used for opaque objects.
   // this should allow for saving pngs with alpha layer.
   color.a = 1.0;
}
