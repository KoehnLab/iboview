
// note: This shader is adapted from freetype-gl's
// demo-distance-field. See:
//    http://code.google.com/p/freetype-gl/

in vec4 vertex_color;
in vec2 tex_coord;
uniform sampler2D FontTexture; // <- it's a signed distance field.

out vec4 color;
const float glyph_center   = 0.50;
       vec3 outline_color  = vec3(0.0,0.0,0.0);
const float outline_center = 0.55;

vec3 fade_base_color(in vec3 color);
void main(void)
{
    vec3 glyph_color = vertex_color.rgb;

    float dist  = texture(FontTexture, tex_coord).r;
    float width = fwidth(dist);
//     float width = 0.10;
    float alpha = smoothstep(glyph_center-width, glyph_center+width, dist);

    float beta = smoothstep(outline_center-width, outline_center+width, dist);
    color = vec4(glyph_color, alpha);

    vec3 rgb = mix(outline_color, color.rgb, beta);
    color = vec4(rgb, max(color.a, beta));
    color.xyz = fade_base_color(color.xyz);

//     color = vec4(glyph_color, alpha);
    if (color.a < 0.005)
      discard;
}
