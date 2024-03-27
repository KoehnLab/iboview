layout(location=0) in vec3 in_Position;
layout(location=1) in vec2 in_Coord2d;
layout(location=2) in vec2 in_Tex;
layout(location=3) in vec4 in_Color;

out vec2 tex_coord;
out vec4 vertex_color;

uniform mat4 mView;
uniform mat4 mProj;

void main()
{
   // transform position and normal
   gl_Position = mProj * (mView * vec4(in_Position, 1.));
   vec4 vView = mView * vec4(in_Position, 1.);
   vView.xy += in_Coord2d;
   gl_Position = mProj * vView;

//    gl_Position = vec4(in_Coord2d.xy + in_Position.xy, 0., 1.);
   vertex_color = in_Color;
   tex_coord = in_Tex;
}
