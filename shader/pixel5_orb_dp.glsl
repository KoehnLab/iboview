layout(location=0) out vec4 color;

// depth of last layer. Only write if current pixel is in front of that.
// (i.e., current z is smaller than last z)
// uniform sampler2DMS Depth1;
uniform sampler2D Depth1;


// layout(early_fragment_tests) in;
// ^- force depth test enable (however, it kinda is the wrong one...).
//    we should only use this for texels which are not discarded.
//    could be re-enabled if using multiple render targets and disabling
//    zwrite (i.e., reversing the roles of Depth0 and Depth1)

vec4 calc_base_color(bool FlipSides);


void main(void)
{
   ivec2 iCoord2d = ivec2(gl_FragCoord.xy);
//    float fDepth0 = texelFetch(Depth1, iCoord2d, gl_SampleID);
   float fDepth0 = texelFetch(Depth1, iCoord2d, 0).r; // has only one channel.
//    float fDepth0 = texelFetch(Depth1, iCoord2d, 0);
// I guess the problem is that the position is wrong?
// gl_FragCoord.z is apparently only the position of the fragment, not of the
// sample.

   if (gl_FragCoord.z < fDepth0) {
      color = calc_base_color(true);
//       color.w = 1.;
   } else {
      discard;
   }

//    //    color = vec4(0.7,0.4,0.4,1.0);
//    float z = gl_FragCoord.z;
//    z = .5 + 100.5 *(fDepth0 - z);
//    color = vec4(z,z,z,1.);
}
