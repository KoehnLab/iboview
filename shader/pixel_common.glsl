// const float LightnessS = 0.2;
// ^- this looks okay for atoms, but the shape perception of orbitals
//    actually does profit from a shinier look.
// const float LightnessS = 0.7;
// const float LightnessD = 0.7;
uniform float ShaderReg0;
uniform float ShaderReg1;
uniform float ShaderReg2;
uniform float ShaderReg3;
uniform float FadeBias;
uniform float FadeWidth;
in vec4 vertex_color;
in vec3 normal;

uniform vec4 DiffuseColor;

vec4 calc_light(in vec3 vNorm, in vec3 vDir, in float Intensity)
{
   vec4 cOut;
   float CosAngle = clamp(dot(vNorm, vDir),0.,1.);
   vec4 cDiffuse = ShaderReg1*pow(CosAngle,ShaderReg0) * DiffuseColor;
   vec4 SpecularColor = vec4(1.,1.,1.,0.);
//    vec4 cSpecular = .5*Lightness*pow(CosAngle,22.) * SpecularColor;
//    vec4 cSpecular = (ShaderReg2*(-.5*pow(CosAngle,16.)+1.2*pow(CosAngle,64.))) * SpecularColor;
   vec4 cSpecular = (ShaderReg2*(ShaderReg3*pow(CosAngle,16.)+1.2*pow(CosAngle,64.))) * SpecularColor;
   // ^- 8). Who can make the shiniest orbitals? MEEEE!!!

   cOut = vertex_color * cDiffuse + cSpecular;
   return Intensity*cOut;
}

vec3 fade_base_color(in vec3 color) {
//    return color;
//    float rz = clamp(8.0 * (gl_FragCoord.z - 0.5) + 0.5, 0., 1.);
// ^- base variant. returns 0.5 at center of coordinate system, if I am not mistaken.
//    float rz = clamp(8.0 * (gl_FragCoord.z - 0.5) + 0.2, 0., 1.);
   float rz = clamp(FadeWidth * (gl_FragCoord.z - 0.5) + FadeBias, 0., 1.);
//    return vec3(rz,rz,rz);
//    return mix(color.xyz, vec3(1.0,1.0,1.0), rz*rz);
   return mix(color.xyz, vec3(1.0,1.0,1.0), rz);
}

vec4 calc_base_color_with_normal(in vec3 in_Normal)
{
   vec3 vNorm = normalize(in_Normal);

   const vec3 vLightDirection0 = vec3(0.50000000,0.50000000,0.70710678);
   const vec3 vLightDirection1 = vec3(-0.43301270,-0.25000000,0.86602540);
   const vec3 vLightDirection2 = vec3(0.43301270,-0.25000000,0.86602540);

   vec4 color = calc_light(vNorm, vLightDirection0, 1.) +
                calc_light(vNorm, vLightDirection1, 0.6) +
                calc_light(vNorm, vLightDirection2, 0.5);
   // ^- that is Mayavi's standard light setup. I like it.
   color[3] /= clamp(abs(vNorm[2]), 0.1, 1.);
//    color.x = pow(color.x, 1.5);
//    color.y = pow(color.y, 1.5);
//    color.z = pow(color.z, 1.5);
   color.xyz = fade_base_color(color.xyz);

   return color;
}

vec4 calc_base_color(bool FlipSides)
{
   if (!FlipSides)
      return calc_base_color_with_normal(normal);
   else {
      vec3 normal1 = normal;
      if (!gl_FrontFacing)
         normal1 *= -1.;
      return calc_base_color_with_normal(normal1);
   }
}
