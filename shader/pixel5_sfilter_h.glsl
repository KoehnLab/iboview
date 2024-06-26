out vec4 color;
in vec4 PixelPos; // generated by vertex5_pass. image coordinates in -1...1 range.

uniform sampler2D SourceImage;

uniform float FilterWeight[5]; // -2, -1, 0, +1, +2
uniform float FilterStepH; // typically: 1./width of source image.
uniform float FilterStepV; // typically: 1./height of source image.

void main(void)
{
   vec2 uv = vec2(.5 * PixelPos.xy + .5);
   color = vec4(0., 0., 0., 0.);
   color += FilterWeight[0] * texture(SourceImage, uv + vec2(-2.0 * FilterStepH, 0.));
   color += FilterWeight[1] * texture(SourceImage, uv + vec2(-1.0 * FilterStepH, 0.));
   color += FilterWeight[2] * texture(SourceImage, uv + vec2(             0., 0.));
   color += FilterWeight[3] * texture(SourceImage, uv + vec2(+1.0 * FilterStepH, 0.));
   color += FilterWeight[4] * texture(SourceImage, uv + vec2(+2.0 * FilterStepH, 0.));
//    color = vec4(FilterStep, FilterWeight[0], FilterWeight[1], 1.);
//    color = vec4(1., .4,.3, 0.7);
}
