layout(location=0) in vec3 in_Position;
layout(location=2) in vec3 in_Normal;
layout(location=3) in vec4 in_Color;
// ^-- note: this layout (0/2/3) coincided with the (implicit?) layout used by
// the legacy functions (non-glVertexAttribPointer-things). No idea what ended
// up on location=1, and anyway this was most likely not portable in any sense.
// Legacy fixed pipeline code has now been removed.

out vec3 normal;
out vec4 vertex_color;

uniform mat4 mView;
uniform mat3 mNorm; // transformation of normal vectors.
uniform mat4 mProj;

void main()
{
   // transform position and normal
   gl_Position = mProj * (mView * vec4(in_Position, 1.));
   normal = vec3(mNorm * in_Normal);

   // just copy color--we do lighting per pixel.
   vertex_color = in_Color;
}

