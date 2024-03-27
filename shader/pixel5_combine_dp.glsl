// uniform sampler2DMS LayerColor; // RGBA output of current layer
// uniform sampler2DMS LayerDepth; // Z-coordinate of current layer.
uniform sampler2D LayerColor; // RGBA output of current layer

out vec4 color;

void main(void)
{
   ivec2 iCoord2d = ivec2(gl_FragCoord.xy);
//    vec4 vLayerColor = texelFetch(LayerColor, iCoord2d, gl_SampleID);
   vec4 vLayerColor = texelFetch(LayerColor, iCoord2d, 0);

   if (abs(vLayerColor.w) < 1e-2) {
      discard;
   } else {
      color = vLayerColor;
   }
}
